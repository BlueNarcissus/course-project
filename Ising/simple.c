#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "ran_uniform.h"

const int P=8;
const int NumberConfig=100000;
const double TempMin=1;
const double TempMax=50;
const double dTemp=0.1;
const int N=P*P;

double D[N+1][N+1];

void randomconfig(int *M0, int *NNC0);
void rmsm();

int main()
{
    int i,j;
    int M, NNC;
    
    FILE* OutputFile=fopen("NormDensState4.txt","w+");
    
    
    // initialize the random number generator with the system time
    InitializeRandomNumberGenerator(time(0));
    
    for(i=0; i<NumberConfig; i++)
    {
        //generate a random configuration and compute M and NNC
        randomconfig(&M, &NNC);
        
        //update the density of state, make a histogram
        D[(M+N)/2][(NNC+2*N)/4]+=1;
    }
    
    for(i=0; i<N+1; i++)
        for(j=0; j<N+1; j++)
        {
            fprintf(OutputFile, "%d %d %lf \n", 2*i-N, 4*j-2*N, D[i][j]/NumberConfig);
        }
    
    //compute the rms magnetization from TempMin to TempMax
    rmsm();
    
    fclose(OutputFile);
    return 0;
}


void randomconfig(int *M0, int *NNC0)
{
    int Matrix[P][P];
    int i,j;
    int M, NNC;
    
    M=0;
    NNC=0;
    
    // generate the random matrix
    for(i=0; i<P; i++)
        for(j=0; j<P; j++)
        {
            if(RandomNumber()<0.5)
                Matrix[i][j]=-1;
            else
                Matrix[i][j]=1;
        }
    
    for(i=0; i<P; i++)
        for(j=0; j<P; j++)
        {
            // compute M and NNC
            M = M + Matrix[i][j];
            NNC = NNC+ Matrix[i][j]*Matrix[(i+1)%P][j] + Matrix[i][j]*Matrix[i][(j+1)%P];
        }
    
    *M0=M;
    *NNC0=NNC;
}


void rmsm()
{
    double T, Mag2, rmsMag, norm;
    int NumberTemp;
    int i,j,k;
    int M, NNC;
    double Mag4Sum, MagError;
    
    FILE* MagFile=fopen("rmsMag4.txt", "w+");
    
    NumberTemp=((TempMax-TempMin)/dTemp)+1;
    
    for(k=0; k<NumberTemp; k++)
    {
        T=TempMin+k*dTemp;
        Mag2=0;
        norm=0;
        Mag4Sum=0;
        for(i=0; i<N+1; i++)
        {
            M=2*i-N;
            for(j=0; j<N+1; j++)
            {
                NNC=4*j-2*N;
                
                Mag2=Mag2+exp(NNC/T)*D[i][j]*M*M;
                Mag4Sum=Mag4Sum+exp(NNC/T)*D[i][j]*M*M*M*M;
                norm=norm+exp(NNC/T)*D[i][j];
            }
        }
        
        //save the normalized rms magnetization
        rmsMag=sqrt(Mag2/norm);
        MagError=sqrt((Mag4Sum/norm-(Mag2/norm)*(Mag2/norm))/(NumberConfig-1));
        MagError=MagError/(2*rmsMag);
        
        fprintf(MagFile, "%lf %lf %lf\n", T, rmsMag/N, MagError/N);
    }
    
    fclose(MagFile);
}

