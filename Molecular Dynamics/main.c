//
//  main.c
//  MolecularDynamics
//
//  Created by Liang Jin on 09/06/15.
//  Copyright (c) 2015 Liang Jin. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ran_uniform.h"

#define PI 3.1415926

const int TailCorrection=1;
const int NumberOfParticles=75;
const double MaxTime=30;
const double TimeStep=0.009;

typedef struct
{
    double x;
    double y;
    double z;
} VECTOR;

VECTOR Positions[NumberOfParticles];
VECTOR Velocity[NumberOfParticles];

double Box;
double CutOff;

double TotalEnergy;
double Potential;
double TotalVirial;
double KineticEnergy;
double Pressure;
double Temperature;

int NumberofCycles;

void EnergyParticle(VECTOR pos,int i,int jb,double *En,double *Vir);
void EnergySystem(void);
void PairForce(VECTOR pos, int i, VECTOR *Force);


int main(void)
{
    int i,j;
    VECTOR ForceOld, ForceNew;
    
    Box=pow(2*NumberOfParticles, 1./3.);
    CutOff=Box/2.;
    
    FILE* InputFile;
    FILE* OutputFile;
    
    InputFile=fopen("newconfiguration.txt","r");
    OutputFile=fopen("data5.txt", "w+");
    
    for(i=0; i<NumberOfParticles; i++)
    {
        fscanf(InputFile, "%lf %lf %lf %lf %lf %lf", &Positions[i].x, &Positions[i].y, &Positions[i].z, &Velocity[i].x, &Velocity[i].y, &Velocity[i].z);
    }
    
    NumberofCycles=MaxTime/TimeStep;
    
    for(j=0; j<NumberofCycles; j++)
    {
        KineticEnergy=0;
        
        for(i=0; i<NumberOfParticles; i++)
        {
            PairForce(Positions[i], i, &ForceOld);
            
            Positions[i].x+=Velocity[i].x*TimeStep+0.5*ForceOld.x*TimeStep*TimeStep;
            Positions[i].y+=Velocity[i].y*TimeStep+0.5*ForceOld.y*TimeStep*TimeStep;
            Positions[i].z+=Velocity[i].z*TimeStep+0.5*ForceOld.z*TimeStep*TimeStep;
            
            // Periodic bourdary condition
            if(Positions[i].x>=Box)
                Positions[i].x-=Box;
            else if (Positions[i].x<0.0)
                Positions[i].x+=Box;
            
            if(Positions[i].y>=Box)
                Positions[i].y-=Box;
            else if(Positions[i].y<0.0)
                Positions[i].y+=Box;
            
            if(Positions[i].z>=Box)
                Positions[i].z-=Box;
            else if(Positions[i].z<0.0)
                Positions[i].z+=Box;
            
            PairForce(Positions[i], i, &ForceNew);
            
            Velocity[i].x+=0.5*(ForceOld.x+ForceNew.x)*TimeStep;
            Velocity[i].y+=0.5*(ForceOld.y+ForceNew.y)*TimeStep;
            Velocity[i].z+=0.5*(ForceOld.z+ForceNew.z)*TimeStep;

            
            KineticEnergy+=0.5*((Velocity[i].x)*(Velocity[i].x)+(Velocity[i].y)*(Velocity[i].y)+(Velocity[i].z)*(Velocity[i].z));
        }
        
        Temperature=2.0*KineticEnergy/(3.0*NumberOfParticles);
        EnergySystem();
        TotalEnergy=Potential+KineticEnergy;
    
        fprintf(OutputFile, "%lf %lf %lf %lf %lf\n", j*TimeStep, TotalEnergy, Pressure, Temperature, Potential);
    }
   
    fclose(InputFile);
    fclose(OutputFile);
    return 0;
}


// Calculate Potential Energy
// Calculate the energy of particle i with particles j=jb, NumberOfparticles
void EnergyParticle(VECTOR pos,int i,int jb,double *En,double *Vir)
{
    double rij,Virij,Enij;
    VECTOR dr;
    int j;
    
    Enij=0.0;
    Virij=0.0;
    
    for(j=jb;j<NumberOfParticles;j++)
    {
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
            
            // Calculate the energy
            if(rij<CutOff)
            {
                Enij+=exp(-rij)/(rij*rij);
                Virij+=exp(-rij)*(rij+2)/(rij*rij);
            }
        }
        
    *En=Enij;
    *Vir=Virij;
    }
}

// Calculates potential energy and pressure
void EnergySystem(void)
{
    double Eni,Viri;
    int i;
    
    Potential=0.0;
    TotalVirial=0.0;
    
    for(i=0;i<NumberOfParticles-1;i++)
    {
        EnergyParticle(Positions[i],i,i+1,&Eni,&Viri);
        Potential+=Eni;
        TotalVirial+=Viri;
    }
    
    Pressure=0.5+TotalVirial/(6*NumberOfParticles);
    
    // Add tail-correction
    if(TailCorrection)
    {
        Potential+=NumberOfParticles*PI*exp(-CutOff);
        Pressure+=PI/6*(CutOff+3)*exp(-CutOff);
    }
}


void PairForce(VECTOR pos, int i, VECTOR *Force)
{
    int j;
    double rij,Ff;
    VECTOR dr,Fij;

    Fij.x=0;
    Fij.y=0;
    Fij.z=0;
    
    for(j=0; j<NumberOfParticles; j++)
    {
        if(j!=i)
        {
            dr.x=pos.x-Positions[j].x;
            dr.y=pos.y-Positions[j].y;
            dr.z=pos.z-Positions[j].z;
            
            // Apply boundary conditions
            dr.x-=Box*rint(dr.x/Box);
            dr.y-=Box*rint(dr.y/Box);
            dr.z-=Box*rint(dr.z/Box);
            
            rij=sqrt((dr.x)*(dr.x)+(dr.y)*(dr.y)+(dr.z)*(dr.z));
            
            // Calculate Force
            if(rij<CutOff)
            {
                Ff=exp(-rij)/(rij*rij)*(2.0/rij+1.0)/rij;
                Fij.x+=Ff*dr.x;
                Fij.y+=Ff*dr.y;
                Fij.z+=Ff*dr.z;
            }
        }
    }
    
    (*Force).x=Fij.x;
    (*Force).y=Fij.y;
    (*Force).z=Fij.z;
}