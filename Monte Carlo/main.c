//  main.c
//  MonteCarlo
//
//  Created by Liang Jin on 20/04/15.
//  Copyright (c) 2015 Liang Jin. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ran_uniform.h"

#define PI 3.1415926

const int TailCorrection=0;
const int NumberOfParticles=50;
const int NumberOfCycles=4000;

typedef struct
{
    double x;
    double y;
    double z;
} VECTOR;

VECTOR Positions[NumberOfParticles];

double Box;
double CutOff;

double TotalEnergy;
double TotalVirial;
double Pressure;

int NumberOfAttempts;
int NumberOfAcceptedMoves;

double MaximumDisplacement;

void Initialization(void);
void EnergyParticle(VECTOR pos,int i,int jb,double *En,double *Vir);
void EnergySystem(void);
void Mcmove(void);


int main()
{
    int k;
    int i;
    
    Box=pow(2*NumberOfParticles, 1.0/3.0);
    CutOff=2.0*Box;
    MaximumDisplacement=Box/2.0;     // change for each step
    
    InitializeRandomNumberGenerator(time(0));
    
    FILE* InputFile=fopen("configuration.txt", "r");
    FILE* OutputFile=fopen("2c_nt.txt","w+");
    
    //Initialization();
    for(i=0;i<NumberOfParticles;i++)
    {
        fscanf(InputFile,"%lf %lf %lf", &Positions[i].x, &Positions[i].y, &Positions[i].z);
    }
    
    for(k=0; k<NumberOfCycles; k++)
    {
        Mcmove();
        EnergySystem();
        fprintf(OutputFile, "%d %lf %lf %lf\n", k, (double)NumberOfAcceptedMoves/NumberOfAttempts,Pressure, TotalEnergy/NumberOfParticles);
    }
    
    fclose(OutputFile);
    return 0;
}

/*
// Initialization: put all the particles randomly in the box
void Initialization(void)
{
    int i;
    for(i=0; i<NumberOfParticles; i++)
    {
        Positions[i].x=RandomNumber()*Box;
        Positions[i].y=RandomNumber()*Box;
        Positions[i].z=RandomNumber()*Box;
    }
}
*/

// Calculate Energy
// Calculate the energy of particle i with particles j=jb, NumberOfparticles
void EnergyParticle(VECTOR pos,int i,int jb,double *En,double *Vir)
{
    double rij,Virij,Enij;
    VECTOR dr;
    VECTOR n;
    int j;
    int k;
    k=rint(CutOff/Box);   // for cut-off equal to 3L/4, 2L
    
    Enij=Virij=0.0;
    
    for(j=jb;j<NumberOfParticles;j++)
    {
        // interaction with the rest of the particles inside the box
        if (j!=i)
        {
            dr.x=pos.x-Positions[j].x;
            dr.y=pos.y-Positions[j].y;
            dr.z=pos.z-Positions[j].z;
            
            // Apply boundary conditions
            dr.x-=Box*rint(dr.x/Box);
            dr.y-=Box*rint(dr.y/Box);
            dr.z-=Box*rint(dr.z/Box);
            
            rij=sqrt((dr.x)*(dr.x)+(dr.y)*(dr.y)+(dr.z)*(dr.z));
            
            // Calculate the interaction with the rest of the particles
            if(rij<CutOff)
            {
                Enij+=exp(-rij)/(rij*rij);
                Virij+=exp(-rij)*(rij+2)/(rij*rij);
            }
        }
        
        // interaction with image particles outside the box
        for(n.x=-k; n.x<=k; (n.x)++)
        {
            for(n.y=-k; n.y<=k; (n.y)++)
            {
                for(n.z=-k; n.z<=k; (n.z)++)
                {
                    if((abs(n.x)+abs(n.y)+abs(n.z))!=0)
                    {
                        //printf("%f %f %f \n", n.x, n.y, n.z);
                        
                        dr.x=pos.x-Positions[j].x;
                        dr.y=pos.y-Positions[j].y;
                        dr.z=pos.z-Positions[j].z;
         
                        // put the particel in the centre
                        dr.x-=Box*rint(dr.x/Box);
                        dr.y-=Box*rint(dr.y/Box);
                        dr.z-=Box*rint(dr.z/Box);
                        
                        // the distance of this particle with the image particles
                        dr.x+=(n.x)*Box;
                        dr.y+=(n.y)*Box;
                        dr.z+=(n.z)*Box;
         
                        rij=sqrt((dr.x)*(dr.x)+(dr.y)*(dr.y)+(dr.z)*(dr.z));
                    
                        if(rij<CutOff)
                        {
                            Enij+=exp(-rij)/(rij*rij);
                            Virij+=exp(-rij)*(rij+2)/(rij*rij);
                        }
                    }
                }
            }
        }

    }
    *En=Enij;
    *Vir=Virij;
}

// Calculates total energy
void EnergySystem(void)
{
    double Eni,Viri;
    int i;
    
    TotalEnergy=0.0;
    TotalVirial=0.0;
    for(i=0;i<NumberOfParticles;i++)
    {
        EnergyParticle(Positions[i],i,i,&Eni,&Viri);
        TotalEnergy+=Eni;
        TotalVirial+=Viri;
    }
    
    Pressure=0.5+TotalVirial/(6*NumberOfParticles);
    
    // Add tail-correction
    if(TailCorrection)
    {
        TotalEnergy+=NumberOfParticles*PI*exp(-CutOff);
        Pressure+=PI/6*(CutOff+3)*exp(-CutOff);
        //printf("%lf %lf\n", NumberOfParticles*PI*exp(-CutOff), PI/6*(CutOff+3)*exp(-CutOff));
    }
    
}


// Monte Carlo Move
// Attempts to displace all the particles in the box
void Mcmove(void)
{
    double EnergyNew,VirialNew,EnergyOld,VirialOld;
    VECTOR NewPosition;
    int i;
    
    for(i=0; i<NumberOfParticles; i++)
    {
    NumberOfAttempts++;
    
    // calculate old energy
    EnergyParticle(Positions[i],i,0,&EnergyOld,&VirialOld);
    
    // give a random displacement
    NewPosition.x=Positions[i].x+(RandomNumber()-0.5)*MaximumDisplacement;
    NewPosition.y=Positions[i].y+(RandomNumber()-0.5)*MaximumDisplacement;
    NewPosition.z=Positions[i].z+(RandomNumber()-0.5)*MaximumDisplacement;
    
    // calculate new energy
    EnergyParticle(NewPosition,i,0,&EnergyNew,&VirialNew);
    
        if(RandomNumber()<exp(-(EnergyNew-EnergyOld)))
        {
         // accept
         NumberOfAcceptedMoves++;
            
         // put particle in simulation box
         if(NewPosition.x<0.0)
         NewPosition.x+=Box;
         else if(NewPosition.x>=Box)
         NewPosition.x-=Box;
         
         if(NewPosition.y<0.0)
         NewPosition.y+=Box;
         else if(NewPosition.y>=Box)
         NewPosition.y-=Box;
         
         if(NewPosition.z<0.0)
         NewPosition.z+=Box;
         else if(NewPosition.z>=Box)
         NewPosition.z-=Box;
        
         // update new position
         Positions[i].x=NewPosition.x;
         Positions[i].y=NewPosition.y;
         Positions[i].z=NewPosition.z;
        }
    }
}
