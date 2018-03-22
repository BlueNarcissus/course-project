//
//  main.c
//  Initialization
//
//  Created by Liang Jin on 09/06/15.
//  Copyright (c) 2015 Liang Jin. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ran_uniform.h"

const int NumberOfParticles=75;

typedef struct
{
    double x;
    double y;
    double z;
} VECTOR;

VECTOR Positions[NumberOfParticles];
VECTOR Velocity[NumberOfParticles];

double Box;
double KineticEnergy;
double Scale;

// Initialization: put all the particles randomly in the box
int main(void)
{
    int i;
    VECTOR Velocity_cm;
    Velocity_cm.x=0.0;
    Velocity_cm.y=0.0;
    Velocity_cm.z=0.0;
    
    InitializeRandomNumberGenerator(time(0));
    FILE* OutputFile=fopen("newconfiguration.txt","w+");
    
    Box=pow(2*NumberOfParticles, 1./3.);
    
    for(i=0; i<NumberOfParticles; i++)
    {
        Positions[i].x=RandomNumber()*Box;
        Positions[i].y=RandomNumber()*Box;
        Positions[i].z=RandomNumber()*Box;
        
        Velocity[i].x=RandomNumber()-0.5;
        Velocity[i].y=RandomNumber()-0.5;
        Velocity[i].z=RandomNumber()-0.5;
        
        Velocity_cm.x+=Velocity[i].x;
        Velocity_cm.y+=Velocity[i].y;
        Velocity_cm.z+=Velocity[i].z;
    }
    
    // Velocity of the Centre of Mass
    Velocity_cm.x=Velocity_cm.x/NumberOfParticles;
    Velocity_cm.y=Velocity_cm.y/NumberOfParticles;
    Velocity_cm.z=Velocity_cm.z/NumberOfParticles;
    
    // Velocity relative to Centre of Mass, the sum of all the velocities are zero
    KineticEnergy=0;
    for(i=0; i<NumberOfParticles; i++)
    {
        Velocity[i].x-=Velocity_cm.x;
        Velocity[i].y-=Velocity_cm.y;
        Velocity[i].z-=Velocity_cm.z;
        
        KineticEnergy+=0.5*((Velocity[i].x)*(Velocity[i].x)+(Velocity[i].y)*(Velocity[i].y)+(Velocity[i].z)*(Velocity[i].z));
    }
    
    Scale=sqrt(NumberOfParticles/KineticEnergy);
    
    KineticEnergy=0;
    Velocity_cm.x=0.0;
    Velocity_cm.y=0.0;
    Velocity_cm.z=0.0;
    
    for(i=0; i<NumberOfParticles; i++)
    {
        Velocity[i].x*=Scale;
        Velocity[i].y*=Scale;
        Velocity[i].z*=Scale;
        
        fprintf(OutputFile, "%lf %lf %lf %lf %lf %lf\n", Positions[i].x, Positions[i].y, Positions[i].z, Velocity[i].x, Velocity[i].y, Velocity[i].z);
    }

    fclose(OutputFile);
    return 0;
}
