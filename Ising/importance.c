#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "ran_uniform.h"

const int P=8;
const int NumberFlip=50000000;
const double T0=2;
const double TempMin=1.5;
const double TempMax=50;
const double dTemp=0.01;
const int N=P*P;

int M, NNC;

int Matrix[P][P];
double D[N+1][N+1];
double N0[N+1][N+1];
double H[N+1];

void randomconfig();
void metropolis();
void Mbar();
void rmsm2();

int main()
{
    int i;
    
    // initialize the random number generator with the system time
    InitializeRandomNumberGenerator(time(0));

    randomconfig();
    
    for(i=0; i<NumberFlip; i++)
    {
        metropolis();
    }
    
    // compute the current rms magnetization
    Mbar();
    
    // compute the density of state and the rms magnetization from TempMin to TempMax
    rmsm2();
    
    return 0;
}


void randomconfig()
{
    int i,j;
    
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
            NNC = NNC + Matrix[i][j]*Matrix[(i+1)%P][j] + Matrix[i][j]*Matrix[i][(j+1)%P];
        }
}


void metropolis()
{
    int ir, jr, im1, jm1, ip1, jp1;
    int Sij;
    double Mn, NNCn;
    double dE;
    
    // choose randomly a spin and flip it
    ir=P*RandomNumber();
    jr=P*RandomNumber();
    
    Sij=Matrix[ir][jr];
    
    // compute dE
    im1=ir-1;
    jm1=jr-1;
    if(im1<0) im1=P-1;
    if(jm1<0) jm1=P-1;
    
    ip1=(ir+1)%P;
    jp1=(jr+1)%P;
    
    Mn=M-Sij-Matrix[ir][jr];
    NNCn=NNC+(-Sij-Matrix[ir][jr])*(Matrix[im1][jr]+Matrix[ir][jm1]+Matrix[ip1][jr]+Matrix[ir][jp1]);;

    dE=-(NNCn-NNC);
    
    // metropolis criterion
    if(RandomNumber()< exp(-dE/T0))
    {
        Matrix[ir][jr]=-Sij;
        M=Mn;
        NNC=NNCn;
    }
    // update the histogram
    N0[(M+N)/2][(NNC+2*N)/4]+=1;
    H[(M+N)/2]+=1;
}

void Mbar()
{
    int i,j;
    int M0, NNC0;
    double Mag2, Mrms, norm;
    double Mag4Sum, MagError;
    
    FILE* OutputFile=fopen("MagT0.txt", "a+");
    FILE* NdosFile=fopen("Ndos_2.txt", "w+");
    
    Mrms=0;
    norm=0;
    Mag4Sum=0;
    
    for(i=0; i<N+1; i++)
    {
        M0=2*i-N;
        for(j=0; j<N+1; j++)
        {
            NNC0=4*j-2*N;
            
            Mag2=Mag2+N0[i][j]*M0*M0;
            Mag4Sum=Mag4Sum+N0[i][j]*M0*M0*M0*M0;
            norm=norm+N0[i][j];
            
            fprintf(NdosFile, "%d %d %lf\n", M0, NNC0, N0[i][j]/NumberFlip);
        }
    }
    
    //save the normalized rms magnetization and the standard error
    Mrms=sqrt(Mag2/norm);
    MagError=sqrt((Mag4Sum/norm-(Mag2/norm)*(Mag2/norm))/(NumberFlip-1));
    MagError=MagError/(2*Mrms);
    
    fprintf(OutputFile, "%lf %lf %lf\n", T0, Mrms/N, MagError/N);
    
    fclose(OutputFile);
    fclose(NdosFile);
}


void rmsm2()
{
    double T, Mag2, Mrms, norm;
    int NumberTemp;
    int i,j,k;
    int M0, NNC0;
    double norm2, Mag4Sum, MagError;
    double Heat, Heat2, HeatCapacity;
    
    FILE* DensFile=fopen("DensState_2.txt", "w+");
    FILE* MagFile=fopen("rmsMag_2.txt", "w+");
    FILE* HeatFile=fopen("heatcap_2.txt", "w+");

    // compute the density of state
    for(i=0; i<N+1; i++)
    {
        M0=2*i-N;
        for(j=0; j<N+1; j++)
        {
            NNC0=4*j-2*N;
            D[i][j]=N0[i][j]*exp(-NNC0/T0);
            
            fprintf(DensFile, "%d %d %lf\n", M0, NNC0, D[i][j]/NumberFlip);
        }
    }
    
    // compute the rms magnetization
    NumberTemp=(TempMax-TempMin)/dTemp+1;
    
    for(k=0; k<NumberTemp; k++)
    {
        T=TempMin+k*dTemp;
        
        Mag2=0;
        norm=0;
        Mag4Sum=0;
        norm2=0;
        
        Heat=0;
        Heat2=0;
        
        for(i=0; i<N+1; i++)
        {
            M0=2*i-N;
            for(j=0; j<N+1; j++)
            {
                NNC0=4*j-2*N;
                
                Mag2+=exp(NNC0/T)*D[i][j]*M0*M0;
                norm+=exp(NNC0/T)*D[i][j];
                Mag4Sum+=exp(NNC0/T)*D[i][j]*M0*M0*M0*M0;
                norm2+=exp(NNC0/T-NNC0/T0)*exp(NNC0/T-NNC0/T0)*N0[i][j];
                
                // calculate heat capacity
                Heat+=exp(NNC0/T)*D[i][j]*NNC0;
                Heat2+=exp(NNC0/T)*D[i][j]*NNC0*NNC0;
            }
        }
        
        // the normalized rms magnetization and the reweighted error
        Mrms=sqrt(Mag2/norm);
        MagError=Mag4Sum/norm-(Mag2/norm)*(Mag2/norm);
        MagError=MagError*norm2/(norm*norm);
        MagError=sqrt(MagError)/(2*Mrms);
        
        HeatCapacity=(Heat2/norm-(Heat*Heat)/(norm*norm))/(T*T);

        fprintf(MagFile, "%lf %lf %lf\n", T, Mrms/N, MagError/N);
        fprintf(HeatFile, "%lf %lf\n", T, HeatCapacity);
    }

    fclose(DensFile);
    fclose(MagFile);
    fclose(HeatFile);
}



