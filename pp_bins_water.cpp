#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>

using namespace std;

int main()
{
    string line,fName, typeName;
    int id,i,j,k,n,m,t,tTime,nAtoms,yStep,currentTimeStep,typ;
    double x,y,z,vx,vy,vz,tau1,tau2,tau3,tau4,tau5,tau6,PE,KE,r,rMag,xJ,yJ,zJ;    
    double TAU1, TAU2, TAU3, p, rho, vAVG, vAVG_kerogen, keKerogen, keMiddle, tempMiddle, tempKerogen;
    double velMiddle, velKerogen, volume2, volume3, rhoN1,rhoN2,rhoN3, rhoKerogen;
    double velMiddle2, velMiddle3, massFlowRate, massFlowRate2, massFlowRate3;
    int nMethaneKerogen, nMethaneMiddle, N1, N2, N3;
    double TAU1b, TAU2b, TAU3b, p2, TAU1c, TAU2c, TAU3c, p3;
    double xLo, xHi, yLo, yHi, zLo, zHi;  
    double Lx, Ly, Lz, volume; 
    double kB = 1.38064852e-23;  
        
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-12;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 ps
    double refVolume = refLength * refLength * refLength;
    double refPressure = 100000;

    //********* INPUT ***********
    
    int nTimeSteps= 5000; // CHANGE, use command line: grep -o 'TIMESTEP' dump_meas.lammpstrj | wc -l

    int skipTimeStep = 3000; // when to start measuring bins -> after steady state

    int tSkip= 1000;
    double deltaT = 0.002;
    double mi;
    
    double binWidth = 2;    // coarse
    double binWidth2 = 0.5; // fine

    // input dump main trj
    ifstream data("dump_main.trj",ios::in);
     
    //**************************
    

    
    int nBins=0, nBins2=0; 

    vector<double> molField;
    vector<double> molFieldO;
    vector<double> molFieldH;
    vector<double> keField;    
    vector<double> keOField;
    vector<double> keHField;
    vector<double> velField;
    vector<double> stressField;
     
    vector<double> molField2;
    vector<double> molFieldO2;
    vector<double> molFieldH2;
    vector<double> keField2;    
    vector<double> keOField2;
    vector<double> keHField2;
    vector<double> velField2;
    vector<double> stressField2;

    int count = 0;
    double nAvSteps = 0.0;
    
    int by;
    int c = 0;
    int dN = floor(nTimeSteps/100);
       
    // loop over all timesteps
    for(t=0;t<nTimeSteps;t++)
    { // loop to extract general sim data e.g current timestep, domain size, etc
        for(n=1;n<10;n++)
        {
            if(n == 2)
            {
                data >> currentTimeStep;             
                c++;    
                if(c >= dN )
                {
                    cout << "currentTimeStep = " << currentTimeStep
                        << "; t = " << t << " [ " 
                        << 100*float(t+1)/float(nTimeSteps) 
                        << "% ]" << endl;  
                        
                    c = 0;
                }
            }
            if(n == 4)
            {
                data >> nAtoms;
//                 cout << "nAtoms = " << nAtoms << endl; 
            }
			if(n == 6)
            {
                data >> xLo >> xHi;
            }

            if(n == 7)
            {
                data >> yLo >> yHi;
            }

            if(n == 8)
            {
                data >> zLo >> zHi;
            }            
            getline(data,line);
        }
        
        if(t == 0) 
        {
            
            // compute domain properties
            Lx = xHi - xLo;
            Ly = yHi - yLo;
            Lz = zHi - zLo;

            nBins = ceil(Lx/binWidth);
            
            int halfnBins = ceil(nBins/2);
            nBins = 2*halfnBins;
            binWidth = Lx/double(nBins);
    


            cout << "binwidth update = " << binWidth << endl;
            cout << "nBins = " << nBins << endl;
            
            molField.resize(nBins, 0.0);
            keField.resize(nBins, 0.0); 
            molFieldO.resize(nBins, 0.0);
            keOField.resize(nBins, 0.0);
            molFieldH.resize(nBins, 0.0);
            keHField.resize(nBins, 0.0);
            velField.resize(nBins, 0.0);
            stressField.resize(nBins, 0.0);
            
            
         
            nBins2 = ceil(Lx/binWidth2); 
            int halfnBins2 = ceil(nBins2/2);
            nBins2 = 2*halfnBins2;
            binWidth2 = Lx/double(nBins2);
    
            cout << "binwidth update = " << binWidth2 << endl;
            cout << "nBins = " << nBins2 << endl;
            
            molField2.resize(nBins2, 0.0);
            keField2.resize(nBins2, 0.0);
            molFieldO2.resize(nBins, 0.0);
            keOField2.resize(nBins, 0.0);
            molFieldH2.resize(nBins, 0.0);
            keHField2.resize(nBins, 0.0);
            velField2.resize(nBins2, 0.0);
            stressField2.resize(nBins2, 0.0);
        } 

        for(n=0;n<nAtoms;n++)
        {
            data>>id>>typ>>x>>y>>z>>vx>>vy>>vz>>tau1>>tau2>>tau3;
            
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

                if(t >= skipTimeStep)
                {
                    if( (x >= xLo) && (x <= xHi) )
                    {
                    
                        by = floor((x-xLo)/binWidth);
                    
                        if(by < 0)
                        {
                            by = 0;
                        }
                        if(by >= nBins)
                        {
                            by = nBins-1;
                        }
                    
                        molField[by] += 1.0;
                        if (typ == 1)
                        {
                            keOField[by] += mi * (vx * vx + vy * vy + vz * vz);
                            molFieldO[by] += 1.0;
                        }

                        if (typ == 2)
                        {
                            keHField[by] += mi * (vx * vx + vy * vy + vz * vz);
                            molFieldH[by] += 1.0;
                        }
                        velField[by] += vx; 
                        stressField[by] += (tau1 + tau2 + tau3);
                        
                        // fine measurements
                       
                        by = floor((x-xLo)/binWidth2);
                    
                        if(by < 0)
                        {
                            by = 0;
                        }
                        if(by >= nBins2)
                        {
                            by = nBins2-1;
                        }
                        
                        molField2[by] += 1.0;
                        if (typ == 1)
                        {
                            keOField2[by] += mi * (vx * vx + vy * vy + vz * vz);
                            molFieldO2[by] += 1.0;
                        }

                        if (typ == 2)
                        {
                            keHField2[by] += mi * (vx * vx + vy * vy + vz * vz);
                            molFieldH2[by] += 1.0;
                        }
                        velField2[by] += vx; 
                        stressField2[by] += (tau1 + tau2 + tau3);
                    }
                }
            }
        }

        if(t >= skipTimeStep)
        {
            nAvSteps += 1.0;
        }       
        
        //****
        getline(data,line);
    }
    
    
    
    // bin measurements   
    
    cout << "bin measurements" << endl;
    
    {
        ofstream binFile("3_bins_meas_liquid.txt",ios::out);
             
        double bin, rho = 0.0, vel =0.0, temp =0.0, binVol, binVol2, pressure;
    
        for(i=0;i<nBins;i++)
        {
            bin = xLo + binWidth*0.5 + binWidth*i;
        
            binVol = binWidth*Ly*Lz*refVolume;

            mi = (15.9994 + 2*1.008) * refMass;

            rho = (molField[i]/3)*mi/(binVol*nAvSteps);
            
            vel = 0.0;
            temp = 0.0;
        
            if(molField[i] > 0.0)
            {
                vel = velField[i]*refVelocity/molField[i];
                temp = (keOField[i] * refVelocity * refVelocity + keHField[i] * refVelocity * refVelocity) / (kB * (2 * (molFieldO[i] + molFieldH[i]) - 3));
            } 
            
            pressure = -(stressField[i]/(3.0*binVol*nAvSteps/refVolume)) *refPressure/1e6;
        
            binFile << bin << '\t'
                    << rho << '\t'
                    << vel << '\t'
                    << temp << '\t'
                    << pressure << '\t'
                    << molField2[i] << '\t'
                    << endl;
                
        }
    }
    

    {
        ofstream binFile("3_bins_meas_fine_liquid.txt",ios::out);
             
        double bin, rho = 0.0, vel =0.0, temp =0.0, binVol, binVol2, pressure;
    
        for(i=0;i<nBins2;i++)
        {
            bin = xLo + binWidth2*0.5 + binWidth2*i;
            
            binVol = binWidth2*Ly*Lz*refVolume;

            mi = (15.9994 + 1.008) * refMass;
            rho = (molField2[i]/3)*mi/(binVol*nAvSteps);
            
            vel = 0.0;
            temp = 0.0;
        
            if(molField2[i] > 0.0)
            {
                vel = velField2[i]*refVelocity/molField2[i];
                temp = (keOField2[i] * refVelocity * refVelocity + keHField2[i] * refVelocity * refVelocity) / (kB * (2 * (molFieldO2[i] + molFieldH2[i]) - 3));
            } 
            
            pressure = -(stressField2[i]/(3.0*binVol*nAvSteps/refVolume)) *refPressure/1e6;
        
            binFile << bin << '\t'
                    << rho << '\t'
                    << vel << '\t'
                    << temp << '\t'
                    << pressure << '\t'
                    << molField2[i] << '\t'
                    << endl;
                
        }
    }


    return 0;
}
