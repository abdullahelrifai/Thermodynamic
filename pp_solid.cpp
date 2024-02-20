#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>

using namespace std;

//SCRIPT TO COMPUTE AVERAGE WALL THERMOMECHANICAL PROPERTIES PER TIMESTEP FOR THE MAIN RUN FILE
//PLOTTING THE OUTPUTTED VALUES VS TIME CAN SHOW THE TIME-EVOLUTION OF THESE PROPERTIES

int main()
{
    string line,fName;
    int id,i,j,k,n,m,t,tTime,nAtoms,yStep,currentTimeStep,typ;

    double x,y,z,vx,vy,vz,tau1,tau2,tau3,tau4,tau5,tau6,PE,KE,r,rMag;
    
    double TAU1Left, TAU2Left, TAU3Left,TAU1Right, TAU2Right, TAU3Right;
    double pistonLeftX, pistonRightX, barrierLeftX, barrierRightX;
    int pistonLeftCount,pistonRightCount,barrierLeftCount,barrierRightCount;
    int nWallRight, nWallLeft;
    double heightLeft, heightRight, volumeLeft, volumeRight, rhoLeft, rhoRight, pLeft, pRight;
    double keTop, keBottom, tempRight, tempLeft;
    double kB = 1.38064852e-23;

    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-12;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 ps
    double refVolume = refLength * refLength * refLength;
    double refPressure = 100000;
    
    //******* INPUT *************
    
    int nTimeSteps=5000; // CHANGE, use command line: grep -o 'TIMESTEP' dump_eq1.lammpstrj | wc -l

    double mi = 196.97 * refMass; // mass of one wall atom

    int tSkip= 1000; // from dump command 
    
    double deltaT = 0.002;

    ifstream data("dump_main.trj",ios::in);
    
    //*************************

    ofstream topFile("3_main_left_vs_time.txt",ios::out);
    ofstream bottomFile("3_main_right_vs_time.txt",ios::out); 

  
    for(t=0;t<nTimeSteps;t++)
    {
        for(n=1;n<10;n++)
        {
            if(n == 4)
            {
                data >> nAtoms;
            }
            
            if(n == 2)
            {
                data >> currentTimeStep;
                
                cout << "currentTimeStep = " << currentTimeStep
                    << "; t = " << t << " [ " 
                    << 100*float(t+1)/float(nTimeSteps) 
                    << "% ]" << endl;    
            }
            
            getline(data,line);
        }
        
        keTop=0.0;
        keBottom=0.0;

        nWallRight = 0;
        nWallLeft = 0;
        
        for(n=0;n<nAtoms;n++)
        {
            data>>id>>typ>>x>>y>>z>>vx>>vy>>vz>>tau1>>tau2>>tau3;
            
            if( (typ == 3) ||  (typ == 4) )
            {
                nWallLeft++;
                keBottom += mi*(vx*vx + vy*vy + vz*vz);
            }
            if( (typ == 5) ||  (typ == 6) )
            {
                nWallRight++;
                keTop += mi*(vx*vx + vy*vy + vz*vz);
            }
        }
        
        tempRight = 0.0;
        
        if(nWallRight > 0)
        {
            tempRight = keTop*refVelocity*refVelocity/(3.0*kB*double(nWallRight));
        }
        
        tempLeft = 0.0;
        if(nWallLeft > 0)
        {
            tempLeft = keBottom*refVelocity*refVelocity/(3.0*kB*double(nWallLeft));
        }
               
        
        getline(data,line);
        
        
        if(t > 0)
        {
            topFile<< deltaT*t*tSkip*refTime <<'\t'
                    << tempLeft <<'\t'
                    << endl;
                
            bottomFile<< deltaT*t*tSkip*refTime <<'\t'
                    << tempRight <<'\t'
                    << endl;      
        }      
                
    }


    return 0;
}
