#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>

//SCRIPT TO COMPUTE AVERAGE LIQUID THERMOMECHANICAL PROPERTIES PER TIMESTEP FOR THE EQUILIBRATION FILE
//PLOTTING THE OUTPUTTED VALUES VS TIME CAN SHOW THE TIME-EVOLUTION OF THESE PROPERTIES, THUS TESTING IF EQUILIBRIUM HAS BEEN REACHED

using namespace std;

int main()
{
    string line, fName;
    int id, n, t, i, j, tTime, nAtoms, yStep, currentTimeStep, typ, index, co, nMols, count;
    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE, r;
    double TAU1Gas, TAU2Gas, TAU3Gas, keGas, velGas, nGas, volumeGas, pGas, rhoGas, tempGas;
    double TAU1Liquid, TAU2Liquid, TAU3Liquid, keLiquid, velLiquid, nLiquid, volumeLiquid, pLiquid, rhoLiquid, tempLiquid, keO, keH, nO, nH;
    double TAU1Solid, TAU2Solid, TAU3Solid, keSolid, velSolid, nSolid, volumeSolid, pSolid, rhoSolid, tempSolid;

    double xLo, xHi, yLo, yHi, zLo, zHi, Lx, Ly, Lz;
    double kB = 1.38064852e-23;

    double refLength = 1e-10;
    double refTime = 1e-12;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 ps
    double refPressure = 100000;
    double refVolume = refLength * refLength * refLength;
    double refMass = 1.6605402e-27;



    //****** MAIN INPUTS *********

    int nTimeSteps = 5000; // CHANGE, use command line: grep -o 'TIMESTEP' dump_equil.lammpstrj | wc -l    //// what is this command and where do you use it

    int tSkip = 1000;
    double deltaT = 0.002;

    double mi;

    double xMin = 80; // we are measuring properties between these limits (xMin - xMax)

    double xMax = 100;


    // import dump trj into "data"
    ifstream data("dump_main.trj", ios::in);

    //***************************

    ofstream liquidFile("3_main_fluid_vs_time.txt", ios::out);

    // loop for all timesteps
    for (t = 0; t < nTimeSteps; t++)
    { 
        for (n = 1; n < 10; n++)
        {
            if (n == 4)
            {
                data >> nAtoms;
            }

            if (n == 2)
            {
                data >> currentTimeStep;

                cout << "currentTimeStep = " << currentTimeStep
                    << "; t = " << t << " [ "
                    << 100 * float(t + 1) / float(nTimeSteps)
                    << "% ]" << endl;
            }

            if (n == 6)
            {
                data >> xLo >> xHi;
            }

            if (n == 7)
            {
                data >> yLo >> yHi;
            }

            if (n == 8)
            {
                data >> zLo >> zHi;
            }

            getline(data, line);
        }


        // compute domain info
        Lx = xMax - xMin;
        Ly = yHi - yLo;
        Lz = zHi - zLo;

        // preinitialise
        TAU1Liquid = 0.0;
        TAU2Liquid = 0.0;
        TAU3Liquid = 0.0;
        keLiquid = 0.0;
        keO = 0.0;
        keH = 0.0;
        nO = 0.0;
        nH = 0.0;
        velLiquid = 0.0;
        nLiquid = 0.0;

        // read atomic data and compute relevant values
        for (n = 0; n < nAtoms; n++) 
        {
            data >> id >> typ >> x >> y >> z >> vx >> vy >> vz >> tau1 >> tau2 >> tau3;

            if (typ == 1 || typ == 2) // fluid   
            {
                if (typ == 1)
                {
                    mi = 15.9994 * refMass; // mass of oxygen atom
                }

                if (typ == 2)
                {
                    mi = 1.008 * refMass; // mass of hydrogen atom
                }

                if ((x >= xMin) && (x <= xMax))
                {
                    
                    TAU1Liquid += tau1;
                    TAU2Liquid += tau2;
                    TAU3Liquid += tau3;
                    if (typ == 1)
                    {
                        keO += mi * (vx * vx + vy * vy + vz * vz);
                        nO += 1.0;
                    }

                    if (typ == 2)
                    {
                        keH += mi * (vx * vx + vy * vy + vz * vz);
                        nH += 1.0;
                    }
                    nLiquid += 1.0;
                    velLiquid += vy;
                }
            }
        }

        getline(data, line);

        //compute relevant properties
        volumeLiquid = Lx * Ly * Lz;

        mi = (15.9994 + 2 * 1.008) * refMass;

        pLiquid = -((TAU1Liquid + TAU2Liquid + TAU3Liquid) / (3.0 * volumeLiquid)) * refPressure / 1e6;
        rhoLiquid = mi * (nLiquid / 3) / (volumeLiquid * refLength * refLength * refLength);

        tempLiquid = 0.0;
        velLiquid = 0.0;

        // average out the properties.
        if (nLiquid > 0.0)
        {
            tempLiquid = (keO * refVelocity * refVelocity + keH * refVelocity * refVelocity) / (kB * (2 * (nO + nH) - 3));
            velLiquid *= refVelocity / nLiquid;
        }

        // output average density, temp, velocity for each timestep.
		
        liquidFile << deltaT * t * tSkip * refTime << '\t'
            << pLiquid << '\t'
            << rhoLiquid << '\t'
            << tempLiquid << '\t'
            << velLiquid << '\t'
            << nLiquid
            << endl;
    }



    return 0;
}
