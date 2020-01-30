#include "parametersAndLibraries.h"
#ifndef __FUNCTIONS_INCLUDED__
#define __FUNCTIONS_INCLUDED__
//Define stationary-code versions of functions if not defined by linear code
class dCell{//dendritic cell
public:
    dCell(){
        x=-1; y=-1; z=-1;
        cogAgRatio=0; probActivation=0; timeAntigenCountLastUpdated=0;
    }

    void set_position(double nx, double ny, double nz){
        x=nx; y=ny; z=nz;
    }
    double get_x() const {return x;}
    double get_y() const {return y;}
    double get_z() const {return z;}

    void set_num_cog_antigen(double t, double cogRatio, double cogPrecision, double *cognateAgProbabilityTable){
        cogAgRatio=cogRatio;
        probActivation=cognateAgProbabilityTable[int(cogRatio/cogPrecision)];
        timeAntigenCountLastUpdated=t;
    }
    double get_cog_antigen_ratio() const { return cogAgRatio; }
    double get_probActivation() const {return probActivation;}

    int update_num_cog_antigen(double time, double timeStep, double antigenDecayRate, double cogPrecision, double *cognateAgProbabilityTable){
        cogAgRatio *= exp(-antigenDecayRate*(time-timeAntigenCountLastUpdated)*timeStep); //Exponential decay
        probActivation = cognateAgProbabilityTable[int(cogAgRatio/cogPrecision)];
        timeAntigenCountLastUpdated=time;
        return probActivation<PROBABILITY_TOLERANCE; //Two birds with one stone
    }
    double get_timeAntigenUpdated() const { return timeAntigenCountLastUpdated; }
    void set_timeAntigenUpdated(double t) {timeAntigenCountLastUpdated=t;}

    int cannot_activate_t_cells() const { return probActivation<PROBABILITY_TOLERANCE; }
private:
    double x,y,z;						//position
    double cogAgRatio, probActivation;//the number of cognate antigen on this DC, and so the probability it will activate a bound Tcell
    double timeAntigenCountLastUpdated; //The last time that the number of cognate antigen on this DC was correct - between then and now, some will have unbound.
};

//Interactions
#define CHECK_CONTACT_WITH_DENDRITES_T(x, y, z, inContact) {																							\
			SET_COORDINATES(x,y,z,coordX,coordY,coordZ);																										\
			nearbyDCs = occupiedPositions[THREED_TO_ONED(coordX,coordY,coordZ)];																				\
			for(auto nbDC: nearbyDCs){																															\
				if(HAS_CONTACTED_DC(tCellNum,nbDC)) continue;/*If the T-cell failed to be activated by this DC before, stop.*/									\
				dendx=dCellList[nbDC].get_x(); dendy=dCellList[nbDC].get_y(); dendz=dCellList[nbDC].get_z();													\
				if( SQMAGNITUDE(x-dendx,y-dendy,z-dendz) <= contactRadius*contactRadius) { dCellNum=nbDC; inContact=1; break;}/*activate T-cell and leave loop.*/\
			}																																					\
		}

#define CHECK_CONTACT_WITH_DENDRITES_DC(x, y, z, inContact) {															\
			SET_COORDINATES(x,y,z,coordX,coordY,coordZ);																		\
			nearbyDCs = occupiedPositions[THREED_TO_ONED(coordX,coordY,coordZ)];													\
			for(auto nbDC: nearbyDCs){																							\
				dendx=dCellList[nbDC].get_x(); dendy=dCellList[nbDC].get_y(); dendz=dCellList[nbDC].get_z();					\
				if( SQMAGNITUDE(x-dendx,y-dendy,z-dendz) <= contactRadius*contactRadius) { dCellNum=nbDC; inContact=1; break;}	\
			}																													\
		}

//DC discrete grid
#define REGENERATE_CELL_POSITIONS_AND_ANTIGEN(newTime,numWavesToGen){																									\
			dCellsPresent=0; for(uint wav=0; wav<numWaves; wav++) {dCellsPresentPerWave[wav]=0; /*nextArrivals[wav]=waveTimeSteps[wav];*/}										\
			for(uint mdci=0; mdci<numDCells*numWavesToGen; mdci++){																												\
				/*Set the number of cognate antigen that this DC has*/																											\
				dCellList[mdci].set_num_cog_antigen(waveTimeSteps[waveNumber], cogAgOnArrival, cognateAgPrecision, cognateAgProbabilityTable);									\
																																												\
				/*Find valid coordinates to place this DC*/																														\
				posFailures=0;																																					\
				/*Keep regenerating coordinates until they are within the sphere and not in contact with another Dendrite*/														\
				while(1==1){ /*More efficient than the better practice of while(inContact==1), because of the need to continue if the random coordinates generated are bad.*/	\
					GENERATE_COORDS; /*Gets random values of x, y and z between -radius and +radius (coordinates from a cube of side = diameter)*/								\
					if(OUTSIDE_SPHERE(x,y,z)) continue;																															\
					inContact=0; CHECK_CONTACT_WITH_DENDRITES_DC(x,y,z,inContact); /*Sets inContact=1 if close to a DC*/														\
					if( (posFailures++)>POSITION_OOB_TOLERANCE ){cerr << "ERROR: Too many attempts to place a DC onto the grid without contacting another. Check parameter choices." << endl; return 1;}\
					else if(inContact==0) break;																																\
				}																																								\
																																												\
				/*If the coordinates are valid, set the DC's position and mark nearby sites on the discrete grid.*/																\
				dCellList[mdci].set_position(x,y,z);																															\
				SET_COORDINATES(x,y,z,coordX,coordY,coordZ);																													\
				uint wav=int(mdci/numDCells);																																	\
				if(newTime>=waveTimeSteps[waveNumber+wav] && timeStepsBetweenDCArrival==0 && (dCellList[mdci].get_probActivation()>0 || noTimeLimit)){ /*This is if we want them all to appear at once*/		\
					PLACE_DC_ON_DISCRETE_GRID(coordX,coordY,coordZ,mdci);																										\
					dCellsPresent++; dCellsPresentPerWave[waveNumber+wav]++;																												\
				}																																								\
			}																																									\
			/*Place T cells on the grid, ensuring they are away from the DCs*/																									\
			posFailures=0;																																						\
			for(auto iter = cellMovementOrder.begin(); iter != cellMovementOrder.end(); iter++){																				\
				tCellNum = *iter;																																				\
				while(1==1){ /*More efficient than the better practice of while(inContact==1), because of the need to continue if the random coordinates generated are bad.*/	\
					GENERATE_COORDS; /*Gets random values of x, y and z between -radius and +radius (coordinates from a cube of side = diameter)*/								\
					if(OUTSIDE_SPHERE(x,y,z)) continue;																															\
					inContact=0; CHECK_CONTACT_WITH_DENDRITES_T(x,y,z,inContact); /*Sets inContact=1 if close to a DC*/															\
					if( (posFailures++)>POSITION_OOB_TOLERANCE ){cerr << "ERROR: Too many attempts to place a T-cell onto the grid without contacting a DC. Check parameter choices." << endl; return 1;}\
					else if(inContact==0) break;																																\
				}																																								\
																																												\
				/*If we have valid coordinates, set them*/																														\
				tCellList[tCellNum].set_position(x, y, z); tCellInitialX[tCellNum]=x; tCellInitialY[tCellNum]=y; tCellInitialZ[tCellNum]=z;										\
				freePathRemaining[tCellNum]=0; tCellList[tCellNum].unset_failed_interaction_ID(); tCellList[tCellNum].reset_num_interactions();									\
			}																																									\
		}

//Definitions for static and linear movement versions of code
class tCell{
public:
    tCell(){//default constructor
        x=-1; y=-1; z=-1;//this must be set later
        vx=0; vy=0; vz=0;//0 velocity
        failedDCContact=-1; interactions=0;
    }

    void set_position(double nx, double ny, double nz){
        x=nx; y=ny; z=nz;
    }
    double get_x() const {return x;}
    double get_y() const {return y;}
    double get_z() const {return z;}
    vector<double> get_pos() const {return {x,y,z};}

    void set_velocity(double nvx, double nvy, double nvz){
        vx=nvx; vy=nvy; vz=nvz;
    }
    double get_vx() const {return vx;}
    double get_vy() const {return vy;}
    double get_vz() const {return vz;}
    vector<double> get_vel() {return {vx,vy,vz};}
    void get_vel_by_ref(vector<double> &v){v[0]=vx; v[1]=vy; v[2]=vz;}

    void set_failed_interaction_ID(int dCellNum) {failedDCContact=dCellNum;}
    void unset_failed_interaction_ID()           {failedDCContact=-1;}
    int get_failed_interaction_ID()              const {return failedDCContact;}
    void reset_num_interactions()                {interactions=0;}
    int increment_num_interactions()             {return (++interactions);} //e.g. if interactions was 0, this returns 1
    int get_num_interactions()                   {return interactions;}
private:
    double x,y,z, vx,vy,vz;
    int failedDCContact=-1; //If a T cell tries to be excited by a DC and fails, this flag prevents it trying every timestep until it has moved away.
    int interactions=0;
};

//Interactions
#define HAS_CONTACTED_DC(tCellNum,dCellNum) (tCellList[tCellNum].get_failed_interaction_ID()==dCellNum)

//DC discrete grid
#define PLACE_DC_ON_DISCRETE_GRID(cx, cy, cz, dcnum) {							\
				for(int p = -1; p<=1; p++){ /*Adjacent cells in x*/							\
					if(cx+p >= numPositions) {break;}										\
					else if(cx+p < 0) {continue;}											\
																							\
					for(int q = -1; q<=1; q++){ /*Adjacent cells in y*/						\
						if(cy+q >= numPositions) {break;}									\
						else if(cy+q < 0) {continue;}										\
																							\
							for(int r = -1; r<=1; r++){ /*Adjacent cells in z*/				\
							if(cz+r >= numPositions) {break;}								\
							else if(cz+r < 0) {continue;}									\
							occupiedPositions[THREED_TO_ONED(cx+p,cy+q,cz+r)].push_back(dcnum);\
						}																	\
					}																		\
				}																			\
			}
#define REMOVE_DC_FROM_DISCRETE_GRID(cx, cy, cz, dcnum){									\
				for(int p = -1; p<=1; p++){ /*Adjacent cells in x*/										\
					if(cx+p >= numPositions) break;														\
					else if(cx+p < 0) continue;															\
																										\
					for(int q = -1; q<=1; q++){ /*Adjacent cells in y*/									\
						if(cy+q >= numPositions) break;													\
						else if(cy+q < 0) continue;														\
																										\
							for(int r = -1; r<=1; r++){ /*Adjacent cells in z*/							\
							if(cz+r >= numPositions) break;												\
							else if(cz+r < 0) continue;													\
							vecIt = &occupiedPositions[THREED_TO_ONED(cx+p,cy+q,cz+r)];		\
							(*vecIt).erase( 															\
								(*vecIt).begin()+distance(												\
									(*vecIt).begin(), find(												\
										(*vecIt).begin(), (*vecIt).end(), dcnum							\
									)/*erases the element 'dcnum' from the vector*/						\
								)																		\
							); 																			\
						}																				\
					}																					\
				}																						\
			}

//Coordinates
#define OUTSIDE_SPHERE(x,y,z) ((x)*(x) + (y)*(y) + (z)*(z) > radius*radius)
#define THREED_TO_ONED(i,j,k) (numPositionsSq*(i) + numPositions*(j) + (k))
#define SET_COORDINATES(x,y,z,coordX,coordY,coordZ) coordX = int((x+radius)/cellSide);\
															coordY = int((y+radius)/cellSide);\
															coordZ = int((z+radius)/cellSide);
#define GENERATE_COORDS {x = genrand_uniform_posneg(engine)*radius; y = genrand_uniform_posneg(engine)*radius; z = genrand_uniform_posneg(engine)*radius;}


//Geometry
#define MAGNITUDE(x,y,z) (sqrt( (x)*(x) + (y)*(y) + (z)*(z) ))
#define SQMAGNITUDE(x,y,z) ((x)*(x) + (y)*(y) + (z)*(z))
#define ROT_MAT(row,col)(rotationMatrix[3*(row)+col])
void arbitraryAxisRotation(double axX, double axY, double axZ, double &vecX, double &vecY, double &vecZ, double angle){
    // double tensorProduct[] = {axX*axX, axX*axY, axX*axZ,
    // 						axY*axX, axY*axY, axY*axZ,
    // 						axZ*axX, axZ*axY, axZ*axZ};
    // double crossProduct[]  = {0, -axZ, axY,
    // 						axZ, 0, -axX,
    // 						-axY, axX, 0};
    double rotationMatrix[] = {	cos(angle) + axX*axX*(1-cos(angle)),	axX*axY*(1-cos(angle))-axZ*sin(angle),	axX*axZ*(1-cos(angle))+axY*sin(angle),
                                   axX*axY*(1-cos(angle))+axZ*sin(angle),	cos(angle) + axY*axY*(1-cos(angle)),	axY*axZ*(1-cos(angle))-axX*sin(angle),
                                   axX*axZ*(1-cos(angle))-axY*sin(angle),	axY*axZ*(1-cos(angle))+axX*sin(angle),	cos(angle)+axZ*axZ*(1-cos(angle))};
    double tempX=vecX, tempY=vecY, tempZ=vecZ;
    vecX = tempX*ROT_MAT(0,0) + tempY*ROT_MAT(0,1) + tempZ*ROT_MAT(0,2);
    vecY = tempX*ROT_MAT(1,0) + tempY*ROT_MAT(1,1) + tempZ*ROT_MAT(1,2);
    vecZ = tempX*ROT_MAT(2,0) + tempY*ROT_MAT(2,1) + tempZ*ROT_MAT(2,2);
}

//Calculation of probabilities

double hyperGeometricFormula(int numAntigenAttachedToDC, int numAntigenInContactArea, int tCellActivationThreshold, double A){
    if(numAntigenAttachedToDC*A<tCellActivationThreshold)	return 0.0;
    if(numAntigenInContactArea<tCellActivationThreshold)	return 0.0;
    if(A==0)	return 0.0;
    if(A==1)	return 1.0;
    boost::math::hypergeometric_distribution<> hypg(A*numAntigenAttachedToDC, numAntigenInContactArea, numAntigenAttachedToDC); //r,n,N
    return 1-cdf(hypg, tCellActivationThreshold-1);
}

//Strings and file IO
#if defined(__linux__)
string get_selfpath() {
    char buff[PATH_MAX]; string sBuff;
    string fname;
    ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
    if (len != -1){
        buff[len] = '\0';
        sBuff=string(buff);
        fname = sBuff.substr(0,sBuff.find_last_of("/"));//The substr is to remove the executable name from the end
        return fname.substr(0,fname.find_last_of("/"));
    }
    else return "\0"; //1 repeat of \0
}
#elif defined (__WIN32) || defined (__CYGWIN__)
string get_selfpath(){
				char buff[PATH_MAX]; string sBuff;
				HMODULE hModule = GetModuleHandleW(NULL);
				ssize_t len= GetModuleFileName(hModule, buff, MAX_PATH);
				if (len != -1){
					buff[len] = '\0';
					sBuff=string(buff);
					return sBuff.substr(0,sBuff.find_last_of("/"));//The substr is to remove the executable name from the end
				}
				else return "\0";
			}
		#else
			#error "Unknown platform"
			return 1
#endif

int fileExists(const char fileName[200]){
    struct stat sb;
    if(stat(fileName, &sb) == 0 && S_ISREG(sb.st_mode)) return 1;
    else return 0;
}

int dirExists(const char directory[150]){
    struct stat sb;
    if(stat(directory, &sb) == 0 && S_ISDIR(sb.st_mode)) return 1;
    else return 0;
}

vector<string> splitString(string inString, char delim, int multiple=1){
    stringstream ss(inString);
    string column; //these columns are split by the delim in the string
    vector<string> outVec;
    int i=0;
    while(getline(ss,column,delim)){
        if(i%multiple==0) outVec.push_back(column);
        i++;
    }
    return outVec;
}

string get_cognate_antigen_ratio_lookup_table_prefix(string configDirectory){
    string configFileName = configDirectory + "/localConfig.txt";
    string prefix = "\0";

    if(!dirExists(configDirectory.c_str())){
        cerr << "Directory " << configDirectory << " not found. Aborting." << endl;
        return "\0";
    }
    else if(!fileExists(configFileName.c_str())){
        cerr << "File " << configFileName << " not found. Create a file containing a line such as " << \
					"COGNATE_ANTIGEN_RATIO_LOOKUP_TABLE_PREFIX /scratch/brown/LNModelPreprocessedProbs/activationProbabilityVersusCognateAntigenRatio" \
					<< endl;
        return "\0";
    }

    else{
        ifstream inFile(configFileName);
        if(!inFile.is_open()){
            cerr << "Could not open file " << configFileName << ". Aborting." << endl;
            return "\0";
        }

        string line; vector<string> cols;
        for(; getline(inFile, line);){

            if(line.npos != line.find("COGNATE_ANTIGEN_RATIO_LOOKUP_TABLE_PREFIX")){
                prefix = line.substr(42); //42 is the length of the tag above, plus a space
                break;
            }

        }

        inFile.close();
        return prefix;
    }
}

string get_cognate_antigen_ratio_lookup_table_prefix(){
    string configDirectory = get_selfpath() + "/config";
    return get_cognate_antigen_ratio_lookup_table_prefix(configDirectory);
}

int createParentDirectories(string path){
    //Get directory names from filename
    vector<string> directoryTree=splitString(path,'/');
    //Edit vector to contain the full path of each directory
    for (uint i=1; i<directoryTree.size()-1; i++) directoryTree[i] = directoryTree[i-1]+"/"+directoryTree[i];
    //Check that each directory exists in turn; create them if not
    for(uint i=0; i<directoryTree.size()-1; i++){ //don't want the last element - that's the filename
        if (directoryTree[i]=="") {continue;}
        else if(!dirExists(directoryTree[i].c_str())){
            try{
                mkdir(directoryTree[i].c_str(),0775); //Using sys/stat. 0755 is due to: rwx rwx r.x = 1+2+4 1+2+4 1+4 = 775
            }catch(const exception& e){
                cerr << "Could not create directory " << directoryTree[i] << endl << e.what() << endl;
                return 1;
            }
        }
    }
    return 0;
}
int createParentDirectories(){
    return createParentDirectories(get_cognate_antigen_ratio_lookup_table_prefix());
}

int generateLookupTable(char outFileName[], int numAntigenAttachedToDC, int numAntigenInContactArea, int tCellActivationThreshold){
    if(createParentDirectories(outFileName)) {return 1;}
    double precision=1E-5; //How many decimals places to store this to
    int numSteps = 1E+5; //NOT relying on division here
    if(numSteps - 1/precision > 1.0/numSteps) {cerr<<"#ERROR: numSteps, " << numSteps << ", not equal to 1/precision, 1/" << precision  << "=" << 1/precision << ". Aborting.\n"; return 1;}
    bool reached_one=false;
    ofstream outFile(outFileName);
    if(!outFile.is_open()){ cerr << "There was a problem opening the file " << outFileName << ". Aborting.\n"; return 1;}
    outFile << "#NUM_ANTIGEN_ON_DC " << numAntigenAttachedToDC << "\n";
    outFile << "#TCELL_ACTIVATION_THRESHOLD " << tCellActivationThreshold << "\n";
    outFile << "#NUM_ANTIGEN_IN_CONTACT_AREA " << numAntigenInContactArea << "\n";
    outFile << "#COGNATE_ANTIGEN_PRECISION " << precision << "\n";
    outFile << "#NUM_VALUES " << numSteps+1 << "\n";
    outFile << "#CognateAntigenRatio ProbabilityOfActivation\n";
    double cogAgR, prob;
    for(int i=0; i<=numSteps; i++){//<= to get 0.0 and 1.0
        cogAgR = i*precision;
        if (!reached_one){
            prob = hyperGeometricFormula(numAntigenAttachedToDC, numAntigenInContactArea, tCellActivationThreshold, cogAgR);
            if (1.0-prob < PROBABILITY_TOLERANCE) reached_one=true;
        }
        else prob=1.0;
        outFile << cogAgR << " " << prob << "\n";
    }

    outFile.close();
    return 0;

}

void lookupTable_printError(string errorReason, char inFileName[]){
    cerr << "Encountered a problem when reading " << inFileName << ". Please investigate. Error: " << errorReason << endl;
    cout << "Simulation aborted due to invalid file " << inFileName << "." << endl;
}

int readLookupTable(char inFileName[], double &precision, int &probabilityTable_size, double **probabilityTable){
    ifstream inFile(inFileName);
    if(!inFile.is_open()){ lookupTable_printError("Failed to open file.",inFileName); return 1;}
    string line; vector<string> cols;
    bool atLeastOneLine=false;

    for(; getline(inFile, line);){
        if(!atLeastOneLine) atLeastOneLine=true;
        cols = splitString(line, ' ');
        if(cols[0].empty()) continue;
        if(cols[0] == "#COGNATE_ANTIGEN_PRECISION") precision = stod(cols[1]);
        else if(cols[0] == "#NUM_VALUES"){
            probabilityTable_size = stoi(cols[1]);
            *probabilityTable = (double*)malloc(probabilityTable_size*sizeof(double));// MALLOC(*probabilityTable, probabilityTable_size, double);
        }
        else if(cols[0] == "#CognateAntigenRatio") break;
    }

    if(!atLeastOneLine){lookupTable_printError("File is empty.",inFileName); return 11;}

    //If we have reached here, then we should be at the actual data table
    getline(inFile,line);
    for(int i=0; i<probabilityTable_size; i++, getline(inFile,line)){
        if(line[0]=='\0') {lookupTable_printError("File ended sooner than expected.",inFileName); return 12;}

        cols=splitString(line, ' ');
        if (cols.size()<2) {lookupTable_printError("Not enough columns in line.",inFileName); return 13;}

        try{
            (*probabilityTable)[i] = stod(splitString(line, ' ')[1]); //The second column, which is the number
        } catch (const invalid_argument&){
            lookupTable_printError("stod: Invalid argument.",inFileName); return 14;
        } catch (const out_of_range&){
            lookupTable_printError("stod: Out of range.",inFileName); return 15;
        }
    }

    inFile.close();
    return 0;
}

string outputDoubleWith1000Seps(double number){
    double divisor=number, i=0; vector<string> places;

    while(int(divisor)!=0){
        divisor /= 1000;  //The number's decimals places are the next 3 target numbers.
        places.push_back( to_string( int(0.5+(divisor-int(divisor))*1000)) ); //Extract the decimal places from the divided number, multiply a thousand, then round up. Convert these 3 digits to a string and store.
        boost::format fmt("%03d"); places.back() = (fmt%places.back()).str(); //Ensure we capture the first leading zero, e.g. 64 -> 064.
        i++;
    }
    //Output result (backwards).
    string returnString="";
    for(auto place=places.rbegin(); place!=places.rend(); place++)  returnString += *place+",";
    return returnString.substr(0,returnString.length()-1);
}

//Vector operations
template<typename T>
vector<int> sort_indexes(const vector<T> &v){
    //Initialise original vector locations
    vector<int> idx(v.size());
    for(size_t i=0; i !=idx.size(); ++i) idx[i] = i;

    //sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2){return v[i1] < v[i2];});
    return idx;
    //source: http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
}
//template <class-type> -- e.g. int, vector
//return-type function-name(parameter-list)

//Rounding
int roundToInt(double input){
    int result = input;
    return result + (input-result>=0.5 ? 1 : 0 );
}
int roundUp(double input){
    if(input-int(input) < TIME_TOLERANCE) return int(input);
    else return int(input+1);
}

//Testing
int criticalASweep(){
    int NValues[] = {DEFAULT_NUM_ANTIGEN_IN_CONTACT_AREA, 2*DEFAULT_NUM_ANTIGEN_IN_CONTACT_AREA, 5*DEFAULT_NUM_ANTIGEN_IN_CONTACT_AREA, 10*DEFAULT_NUM_ANTIGEN_IN_CONTACT_AREA, int(0.1*DEFAULT_NUM_ANTIGEN_IN_CONTACT_AREA)};
    int TValues[] = {DEFAULT_TCELL_ACTIVATION_THRESHOLD, 10*DEFAULT_TCELL_ACTIVATION_THRESHOLD, int(0.1*DEFAULT_TCELL_ACTIVATION_THRESHOLD)};
    vector<vector<double> > criticalAValues(5, vector<double>(3)); //5x3 vector
    double precision=1E-3; //How many decimals places to store this to
    int numSteps = 1E+3; //NOT relying on division here
    int N, T; long int i=0;
    double cogAgR, prob=0;

    for(int n=0; n<5; n++){
        N = NValues[n];
        for(int t=0; t<3; t++){
            cerr << n << "," << t << "\t";
            T = TValues[t]; prob=0;
            criticalAValues[n][t]=1; //In case the value is not set in this loop, it will be "1"

            for(i=0, precision=1E-3, numSteps=1E+3; i<numSteps; i++){
                // if(i%(numSteps/100)==0) cerr << precision << ":" << 100.0*i/numSteps << "% ";
                cogAgR = i*precision;
                prob = hyperGeometricFormula(2000.0/7.8*N, N, T, cogAgR);

                if(prob>=PROBABILITY_TOLERANCE){
                    i--; precision/=10; numSteps*=10; i*=10;
                    if(numSteps>1E+9){
                        cerr << i << endl;
                        criticalAValues[n][t] = i*precision;
                        break;
                    }
                }
            }
        }
    }
    cout << "\t"; for(int n=0; n<5; n++) cout << "N=" << NValues[n] << "\t";
    cout << endl;
    for(int t=0; t<3; t++){
        cout << "T=" << TValues[t] << "\t";
        for(int n=0; n<5; n++) cout << criticalAValues[n][t] << "\t";
        cout << endl;
    }

    return 0;
}
#endif
//__FUNCTIONS_INCLUDED__