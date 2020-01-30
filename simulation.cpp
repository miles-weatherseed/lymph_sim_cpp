#include "parametersAndLibraries.h"
#include "classesAndFunctions.h"
#include "csv.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>

//Prototype (with variable names intact for sanity)
int simulation(long unsigned int numTimeSteps, double timeStep, uint
numTimeMeasurements, vector<double> waveTimes, double radius,  double contactRadius, double tCellAgSpecFreq, uint
               numTCells, uint numDCells, int numAntigenInContactArea,  int numAntigenOnDC, int tCellActivationThreshold, int
               cognateAgProbabilityTable_size, double cognateAgPrecision, double *cognateAgProbabilityTable, double cogAgInDermis,
               double cogAgOnArrival,  unsigned int firstDCArrival, unsigned int DCArrivalDuration, double antigenDecayRate, double
               tGammaShape, double tGammaScale, double freePathMean, double
               freePathStDev,  double dCellVelocity, int seed, int numRepeats, bool noTimeLimit, int *numActivationsPerRepeat, int print);

//Entry point for c++
int loop(double contactRadius=DEFAULT_CONTACT_RADIUS, uint numDCells=DEFAULT_DENDNUM, double freePathMean=DEFAULT_T_FREE_PATH_MEAN,
        double freePathStDev=DEFAULT_T_FREE_PATH_STDEV, int tCellActivationThreshold=DEFAULT_TCELL_ACTIVATION_THRESHOLD,
        double firstDCArrival=DEFAULT_FIRST_DC_ARRIVAL, double DCArrivalDuration=DEFAULT_DC_ARRIVAL_DURATION,
        double radius=DEFAULT_RADIUS, double tVelocityMean=DEFAULT_T_VELOCITY_MEAN, double cogAgInDermis=DEFAULT_COGNATE_RATIO_MEAN_AT_DERMIS,
        double antigenDecayRate=DEFAULT_ANTIGEN_DECAY_RATE, int numAntigenInContactArea=DEFAULT_NUM_ANTIGEN_IN_CONTACT_AREA){
    //Initialise params
    double timeStep = DEFAULT_TIMESTEP;
    long unsigned int numTimeSteps = NUM_TIMESTEPS_2DAY(timeStep);
    uint numTimeMeasurements=DEFAULT_NUM_TIME_MEASUREMENTS;
    double tCellAgSpecFreq = DEFAULT_AGSPEC_FREQ; //n.b. these are antigen-specific T-cells, rho*phi, where rho is the density of all T-cells and phi is the number specific to the antigen of interest
    uint numTCells = TOTAL_TNUM * tCellAgSpecFreq;
    int numAntigenOnDC = numAntigenInContactArea*2000.0/7.8;
    double cogAgOnArrival = DEFAULT_COGNATE_RATIO_MEAN_AT_DERMIS;
    // double cognateAntigenRatioStDev = DEFAULT_COGNATE_RATIO_STDEV;
    unsigned int firstDCArrival = firstDCArrival;
    unsigned int DCArrivalDuration = DCArrivalDuration;
    double tVelocityMean=DEFAULT_T_VELOCITY_MEAN;
    double tVelocityStDev=DEFAULT_T_VELOCITY_STDEV;
    double tGammaShape=T_GAMMA_SHAPE;
    double tGammaScale= tVelocityMean/T_GAMMA_SHAPE;
    double freePathMean=freePathMean;
    double freePathStDev=freePathStDev;
    srand(time(NULL));
    double seed = rand();
    int numRepeats = DEFAULT_NUM_REPEATS;
    double antigenDecayRate = antigenDecayRate;
    string timeStepOpt="unset";
    bool TCellsSet=false, initialAntigenSet=false, noTimeLimit=false;
    vector<double> waveTimes = {0.0};

    // set amount of cogAg on arrival in LN
    cogAgOnArrival=cogAgInDermis*exp(-antigenDecayRate*firstDCArrival);

    //Find lookup table for the probability of activation given ratio of cognate antigen; make one if it is not present
    int cognateAgProbabilityTable_size; double cognateAgPrecision, *cognateAgProbabilityTable;
    string cognate_antigen_ratio_lookup_table_prefix = get_cognate_antigen_ratio_lookup_table_prefix();
    if (cognate_antigen_ratio_lookup_table_prefix=="\0") return 1; //Invalid file or directory -> exit
    char CALTFilename[400]; snprintf(CALTFilename, sizeof(CALTFilename), "%s_%d_%d.dat", cognate_antigen_ratio_lookup_table_prefix.c_str(), numAntigenInContactArea, tCellActivationThreshold);
    if(fileExists(CALTFilename)){
        int retVal;
        if((retVal=readLookupTable(CALTFilename, cognateAgPrecision, cognateAgProbabilityTable_size, &cognateAgProbabilityTable))) return retVal; //returns 1 if an error
    }
    else{
        int retVal;
        if((retVal=generateLookupTable(CALTFilename, numAntigenOnDC, numAntigenInContactArea, tCellActivationThreshold))) return retVal;
        if((retVal=readLookupTable(CALTFilename, cognateAgPrecision, cognateAgProbabilityTable_size, &cognateAgProbabilityTable))) return retVal; //A waste of time, but the time involved is miniscule so I care not
    }

    int numActivationsPerRepeat[numRepeats];

    int returnValue = simulation(numTimeSteps, timeStep, numTimeMeasurements, waveTimes, radius,
                                 contactRadius, tCellAgSpecFreq, numTCells,  numDCells, numAntigenInContactArea, numAntigenOnDC,
                                 tCellActivationThreshold, cognateAgProbabilityTable_size, cognateAgPrecision, cognateAgProbabilityTable,
                                 cogAgInDermis, cogAgOnArrival, firstDCArrival, DCArrivalDuration, antigenDecayRate, tGammaShape, tGammaScale,
                                 freePathMean, freePathStDev, 0, seed, numRepeats, noTimeLimit,  numActivationsPerRepeat, 1);

    //Free all memory we took
    free(cognateAgProbabilityTable);

    return returnValue;
}

int main(){
    io::CSVReader<12> in("param_values.csv");
    in.read_header(0, "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l");
    double a, b, c, d, e, f, g, h, i, j, k, l;
    // double b3, uint D, double F, double f, int n, double P, double p, double R3, double V, double Ad, double k, int N
    while(in.read_row(a, b, c, d, e, f, g, h, i, j, k, l)){
        int val = loop(pow(a, 0.33333), uint(b), c, d, int(e), f, g, pow(h, 0.3333333), i, j, k, int(l));
        cout << val << endl;
    }
}

#include "simulation.h"