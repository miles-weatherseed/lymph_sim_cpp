#include "parametersAndLibraries.h"
#include "classesAndFunctions.h"

int simulation(long unsigned int numTimeSteps, double timeStep, uint
numTimeMeasurements, vector<double> waveTimes, double radius, double contactRadius, double tCellAgSpecFreq, uint
               numTCells, uint numDCells, int numAntigenInContactArea, int numAntigenOnDC, int tCellActivationThreshold, int
               cognateAgProbabilityTable_size, double cognateAgPrecision, double *cognateAgProbabilityTable, double cogAgInDermis,
               double cogAgOnArrival, unsigned int firstDCArrival, unsigned int DCArrivalDuration, double antigenDecayRate, double
               tGammaShape, double tGammaScale, double freePathMean, double freePathStDev,
               double dCellVelocity, int seed, int numRepeats, bool noTimeLimit, int
               *numActivationsPerRepeat, int print=1){

    //More parameters used to hold data later
    double x, y, z, dx, dy, dz, newx, newy, newz, dendx, dendy, dendz, vmag, theta, phi, magsq, dist, minimumCognateAg=0, maxTime=numTimeSteps*timeStep;
    int tCellNum=0, dCellNum=0, dCellsPresent=0, coordX, coordY, coordZ, dcoordX, dcoordY, dcoordZ, inContact=0, numPositions, numPositionsSq;
    uint activationTimeLimit=numTimeSteps, maxSimultaneousDCs=numDCells, numWaves=1, waveNumber=0, tot_numActivated=0, posFailures=0, this_numActivated=0, this_numUniqueInteractions=0, this_numInteractions=0;
    vector<uint> waveTimeSteps(1); waveTimeSteps[0]=int(waveTimes[0]);
    vector<uint> dCellsPresentPerWave(1); dCellsPresentPerWave[0]=0;
    unsigned long int t=0;
    vector<double> radialVec, parallel, reflectedVec, position, velocity(3);
    vector<int> nearbyDCs; vector<int> *vecIt;


    //Find the cognate ag ratio corresponding effectively zero probability of activation
    for(int i=0; i<cognateAgProbabilityTable_size; i++) if(cognateAgProbabilityTable[i]<PROBABILITY_TOLERANCE) minimumCognateAg=cognateAgPrecision*i;
    //Calculate time that the LN cannot activate any T cells anymore.
    if(cogAgOnArrival>minimumCognateAg && minimumCognateAg!=0 && antigenDecayRate!=0 && !noTimeLimit) activationTimeLimit = int(log(cogAgOnArrival/minimumCognateAg)/(antigenDecayRate*timeStep));
    else if(cogAgOnArrival<minimumCognateAg && !noTimeLimit) activationTimeLimit = 0;
    else                                                     activationTimeLimit = numTimeSteps;

    //If multiple doses are being administered, check if any of them will overlap. If so, the memory required to place all the DCs at once will be increased.
    maxSimultaneousDCs=numDCells;

    //Create containers for our T cells and DCs
    vector<tCell> tCellList(numTCells);
    vector<dCell> dCellList(maxSimultaneousDCs);

    //Create a discrete grid to reduce the search space for T cells looking for DCs. However, make this grid iteratively more coarse until the memory usage is less than 1GB.
    long int occupiedPositionsArraySize=0; double contractGridReductionFactor=2.0;//actually starts at 1.0 due to division
    do{
        //Reduce memory usage of grid
        if(contractGridReductionFactor < MINIMUM_CONTACT_GRID_MEMORY_REDUCTION_FACTOR){
            cerr << "Could not produce a discrete position grid of small enough footprint. Check parameters." << endl
                 << "\nnumPositions: " << radius << "\ncontactRadius: " << contactRadius << "\nFull memory usage: " <<  10*(radius / contactRadius+1)*sizeof(int)*1E-6 << "MB"
                 << "\nMinmimum reduction factor: " << MINIMUM_CONTACT_GRID_MEMORY_REDUCTION_FACTOR << "\nMinmum memory usage: " << 10*sizeof(int)*occupiedPositionsArraySize*1E-6 << "MB" << endl;//The 10 estimates how many DCs might take a spot, and the overhead
            return 1;
        }
        else contractGridReductionFactor/=2;

        //Set number of positions and calculate memory usage.
        numPositions = 2*contractGridReductionFactor*radius / contactRadius+1; //e.g. radius = 4, contact radius = 2 -> there are 2*4/2 +1 = 5 positions per dimension, at x=0,2,4,-2,-4
        numPositionsSq = numPositions*numPositions;
        occupiedPositionsArraySize = numPositions*numPositions*numPositions;
    } while(occupiedPositionsArraySize*10*sizeof(int)>2E9);
    double cellSide=contactRadius/contractGridReductionFactor;
    vector<vector<int> > occupiedPositions(occupiedPositionsArraySize);

    //More containers for T cells.
    vector<int> cellMovementOrder(numTCells);
    vector<double> freePathRemaining(numTCells);
    if (numDCells==1) DCArrivalDuration=0;
    double timeStepsBetweenDCArrival = (numDCells>1) ? 1.0*(DCArrivalDuration)/timeStep/(numDCells-1) : 0;//e.g. with a time gap of 6 hours, 11 DCs will arrive at 0, 0.6, 1.2, ... 5.4, 6.0.
    vector<double> tCellInitialX(numTCells),tCellInitialY(numTCells),tCellInitialZ(numTCells);

    //Start random num generator
    mt19937_64 engine(seed); //Random number generator (Mersenne Twister) initialised with seed
    uniform_int_distribution<int> genrand_int(0,1); //generates 1 or 0
    uniform_real_distribution<double> genrand_uniform_posneg(-1,1); //generates doubles between -1 and 1
    uniform_real_distribution<double> genrand_uniform_pos(0,1);
    gamma_distribution<double> genrand_t_velocity(tGammaShape,tGammaScale);
    normal_distribution<double> genrand_t_freepath(freePathMean, freePathStDev);

        for(int repeat=0; repeat<numRepeats; repeat++) {

            //Reinitialise number of activated t-cells and positions of DCs
            this_numActivated = 0;
            this_numUniqueInteractions = 0;
            this_numInteractions = 0;
            waveNumber = 0;
            if (repeat > 0)
                for (int i = 0; i < occupiedPositionsArraySize; i++)
                    occupiedPositions[i].clear(); // = emptyVector; //clear() is marginally faster than the equality. This loop still makes it slow, however.

            if (activationTimeLimit > 0) {
                //Initialise a vector of integers that we will use to randomise movement of the T-cells, and a vector that indicates how long a T cell will move before it "collides"
                cellMovementOrder.resize(numTCells);
                for (uint i = 0; i < numTCells; i++)
                    cellMovementOrder[i] = i; //No need to randomise for first use, since the cells themselves were randomly generated

                //Place cells at random in the sphere and give them a random number of antigen
                REGENERATE_CELL_POSITIONS_AND_ANTIGEN(0, int(maxSimultaneousDCs/numDCells)); //also sets dCellsPresent=0
            } else { //If activation is impossible at the first time step, set the cellMovementOrder vector to 0 such that we record data and exit immediately
                cellMovementOrder.resize(0);
            }

            //Move cells until end of time
            for (t = 0; t < numTimeSteps; t++) {
                //Break if we have exhausted all of the t-cells or if the last wave of DCs are out of cognate antigen
                if (cellMovementOrder.size() == 0 || t > activationTimeLimit) {
                    break;
                }

                //Spawn new DCs
                while ((t >= 0 && t * timeStep <= DCArrivalDuration && DCArrivalDuration != 0 &&
                        t >= dCellsPresentPerWave[0] * timeStepsBetweenDCArrival) ||
                       (timeStepsBetweenDCArrival == 0 && dCellsPresentPerWave[0] < numDCells)) {
                    dCellNum = dCellsPresentPerWave[0];
                    dCellList[dCellNum].set_timeAntigenUpdated(0.0);
                    x = dCellList[dCellNum].get_x();
                    y = dCellList[dCellNum].get_y();
                    z = dCellList[dCellNum].get_z(); //Note that the non-vector version of dCellNum is used here, because we have already ensured that each cell should get a unique position
                    SET_COORDINATES(x, y, z, coordX, coordY, coordZ);
                    if (dCellList[dCellNum].get_probActivation() > 0 || noTimeLimit) {
                        PLACE_DC_ON_DISCRETE_GRID(coordX, coordY, coordZ, dCellNum);
                    }
                    dCellsPresent++;
                    dCellsPresentPerWave[0]++;
                }

                for (auto iter = cellMovementOrder.begin(); iter != cellMovementOrder.end(); iter++) {
                    tCellNum = *iter;
                    //Get current position and record how far it has gone from the beginning
                    x = tCellList[tCellNum].get_x();
                    y = tCellList[tCellNum].get_y();
                    z = tCellList[tCellNum].get_z(); //position = tCellList[tCellNum].get_pos(); x=position[0]; y=position[1]; z=position[2];
                    if (t % (numTimeSteps / numTimeMeasurements) == 0) {
                        dist = MAGNITUDE(x - tCellInitialX[tCellNum], y - tCellInitialY[tCellNum],
                                         z - tCellInitialZ[tCellNum]);
                    }

                    //Move the cell: get current velocity & add it to position
                    velocity = tCellList[tCellNum].get_vel();
                    dx = velocity[0] * timeStep;
                    dy = velocity[1] * timeStep;
                    dz = velocity[2] * timeStep;
                    newx = x + dx;
                    newy = y + dy;
                    newz = z + dz;
                    freePathRemaining[tCellNum] -= MAGNITUDE(dx, dy, dz);

                    //Check if we have tried to left the sphere or the T cell has reached the end of its mean free path
                    if (OUTSIDE_SPHERE(newx, newy, newz)) {
                        //Regenerate velocity to move the T cell away from the sphere surface.
                        vmag = genrand_t_velocity(engine);
                        theta = acos(genrand_uniform_pos(engine));
                        phi = 6.283185307 * genrand_uniform_pos(engine); //random velocity; theta[0,pi/2] and phi[0,2pi]
                        magsq = MAGNITUDE(newx, newy, newz); //Get magnitude of position vector to normalise
                        x /= magsq;
                        y /= magsq;
                        z /= magsq; //Set as a unit vector
                        velocity[0] = -vmag * x;
                        velocity[1] = -vmag * y;
                        velocity[2] = -vmag *
                                      z; //Set velocity direction opposite to the position vector (i.e. pointing to centre of the sphere)
                        arbitraryAxisRotation(z, y, -x, velocity[0], velocity[1], velocity[2],
                                              theta); //Choosing the radial axis as the z axis, this performs a rotation by theta about the x axis.
                        arbitraryAxisRotation(x, y, z, velocity[0], velocity[1], velocity[2],
                                              phi); //Using the same axes as before, this performs a rotation by phi about the z axis.
                        tCellList[tCellNum].set_velocity(velocity[0], velocity[1], velocity[2]);
                        freePathRemaining[tCellNum] = genrand_t_freepath(
                                engine); //Set how long until the particle stops and must choose a new velocity

                        //Finally, reverse the previous step we tried to make and take a new one
                        newx += velocity[0] * timeStep - dx;
                        newy += velocity[1] * timeStep - dy;
                        newz += velocity[2] * timeStep - dz;
                    } else if (freePathRemaining[tCellNum] <=
                               0) {//Check if the cell has exceeded the free path of its previous velocity vector and needs a new one generated
                        vmag = genrand_t_velocity(engine);
                        theta = acos(2 * genrand_uniform_pos(engine) - 1);
                        phi = 6.283185307 * genrand_uniform_pos(engine); //theta between 0 and pi, phi between 0 and 2pi
                        tCellList[tCellNum].set_velocity(vmag * sin(theta) * cos(phi), vmag * sin(theta) * sin(phi),
                                                         vmag * cos(theta)); //Set velocity
                        freePathRemaining[tCellNum] = genrand_t_freepath(
                                engine); //Set how long until the particle stops and must choose a new velocity
                    }

                    //Check for DCs nearby to activate T-cell and either update position or remove from simulation
                    inContact = 0;
                    CHECK_CONTACT_WITH_DENDRITES_T(newx, newy, newz, inContact);
                    if (inContact ==
                        1) {//If this T cell encountered some dendrites, check to see if the interaction was successful
                        //First, update the number of antigen on this DC, as some will have unbound since we last saw the DC.
                        dCellList[dCellNum].update_num_cog_antigen(t, timeStep, antigenDecayRate, cognateAgPrecision,
                                                                   cognateAgProbabilityTable);
                        SET_COORDINATES(dendx, dendy, dendz, dcoordX, dcoordY, dcoordZ);

                        this_numInteractions++;//NOTE: This will be an underestimate if we remove DCs from the discrete grid when they cannot activate T-cells. Note also that T-cells with a "failed interaction ID" with a given DC do not set inContact=1.
                        if (1 ==
                            tCellList[tCellNum].increment_num_interactions()) { this_numUniqueInteractions++; } //increment_num_interactions returns the new value after incrementation

                        if (dCellList[dCellNum].cannot_activate_t_cells() && !noTimeLimit) {
                            REMOVE_DC_FROM_DISCRETE_GRID(dcoordX, dcoordY, dcoordZ,
                                                         dCellNum); //If there is no longer enough cognate Ag, there is little point simulating the DC anymore...
                        }

                            //Now, see if the T cell has been activated
                        else if (genrand_uniform_pos(engine) < dCellList[dCellNum].get_probActivation()) {
                            //Activated: we no longer want to track this cell. Erase it and pull the iterator back to not accidentally skip the next T cell
                            cellMovementOrder.erase(
                                    cellMovementOrder.begin() + distance(cellMovementOrder.begin(), iter));
                            iter--;
                            //Finally, record date for analysis
                            this_numActivated++;
                        } else { //Activation failed - mark this T cell-DC pair so that they won't repeatedly try to interact.
                            tCellList[tCellNum].set_failed_interaction_ID(dCellNum);
                            //Regenerate velocity to move the T cell away from the DC after interaction.
                            vmag = genrand_t_velocity(engine);
                            theta = acos(genrand_uniform_pos(engine));
                            phi = 6.283185307 *
                                  genrand_uniform_pos(engine); //random velocity; theta[0,pi/2] and phi[0,2pi]
                            x = newx - dendx;
                            y = newy - dendy;
                            z = newz - dendz; //Defines the position vector from DC to T cell.
                            magsq = MAGNITUDE(x, y, z); //Get magnitude (not squared) to normalise
                            x /= magsq;
                            y /= magsq;
                            z /= magsq; //Set as a unit vector
                            velocity[0] = vmag * x;
                            velocity[1] = vmag * y;
                            velocity[2] = vmag * z; //Magnitude of velocity, currently parallel to position vector
                            arbitraryAxisRotation(z, y, -x, velocity[0], velocity[1], velocity[2],
                                                  theta); //Setting position vector as z axis, this performs a rotation by theta about the x axis.
                            arbitraryAxisRotation(x, y, z, velocity[0], velocity[1], velocity[2],
                                                  phi); //Using the same axes as before, this performs a rotation by phi about the z axis.
                            tCellList[tCellNum].set_velocity(velocity[0], velocity[1], velocity[2]);
                            freePathRemaining[tCellNum] = genrand_t_freepath(
                                    engine); //Set how long until the particle stops and must choose a new velocity
                        }
                    }

                    //Update cell position (including for ones just removed)
                    tCellList[tCellNum].set_position(newx, newy, newz);

                }//End of T-cell loop

            }//End of time loop

            //Record data for final time
            numActivationsPerRepeat[repeat] = this_numActivated;
            if (this_numActivated > 0) {
                tot_numActivated += 1;
            }
        }//End of repeat loop

    return tot_numActivated;
}