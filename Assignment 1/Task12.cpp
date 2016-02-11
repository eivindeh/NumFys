#include <fstream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <random>
#include "functions.h"
//#define MAXPNT  1000000
using namespace std;

double k_bT      =  26;
double delta_U   =  80000;
double alpha     =  0.2;
double L         =  20;
double tau       =  70;
double dt        =  0.001;
int counter      =  0;
double r         =  12;
double eta       =  1;
const long int N =  800000;
const int M      =  200;

double gamma_i   =  6*3.1415*r*eta;
double D_hat     =  k_bT/delta_U;
double omega     =  delta_U/(gamma_i*L*L);

double x_prev[3][M];
double x_next[3][M];
double chi;


void U_r(double x_hat,double t_hat,double * U)

{   double x_mod=abs(x_hat-floor(x_hat));
    if (x_mod > 0 && x_mod < alpha)
        *U = x_mod*k_bT/(alpha*D_hat);
    else
        *U = (1-x_mod)/(1-alpha)*k_bT/D_hat;
    
    if (t_hat > 0 && t_hat < 3*tau*omega/4)
        *U = 0;
      
}

void F_r(double x_hat,double t_hat,double * F)

{   
    double x_mod=abs(x_hat-floor(x_hat));
    double t_mod = abs(t_hat-floor(t_hat/(tau))*tau);

    if (t_mod >= 0 && t_mod < 3*tau/4)
        *F = 0;

    else if (x_mod >= 0 && x_mod < alpha)
        *F = -1/(alpha);
    else
        *F = 1/(1-alpha);
   
}

double  getRandom(std::default_random_engine generator, std::normal_distribution<double> distribution)
{
  double number = distribution(generator);
  cout<<number;
  return number;
}

double step(double x_n,double t_n, double chi)
{
    double F=0;
    F_r(x_n,t_n,&F);

    //double chi = getRandom();

    return  x_n + F*dt + sqrt(2*D_hat*dt)*chi;
    
}


void print2file(double x[][M],FILE *fptr2)
{
    for (int j = 0; j < M; j++)
    {   
        fprintf(fptr2,"%15.6f%15.6f\n",x[0][j],x[1][j]);
    }
}

double getAvgVel(double x[][M])
{
    double avg;
    for( int j = 0; j<M ; j++){
        avg = avg + x[0][j];
    }
    return avg/x[1][0]/M;
}


int main()

{
    FILE *fptr;
    fptr = fopen("velData.txt","w");
    
    FILE *fptr2;
    fptr2 = fopen("simOut.txt","w");

    cout<<omega<<endl;
    double LHS = dt/alpha+4*sqrt(2*D_hat*dt);
    if (LHS>alpha/10){
        cout<<"Error: Too large timestep"<<endl;
    }

    srand (time(NULL));

    std::default_random_engine        generator;
    std::normal_distribution<double>  distribution(0,1);

for (r=12 ; r <= 36; r=r+24){
    gamma_i   = 6*3.1415*r*eta;
    tau=12*70/r;
    for (int j = 0; j < M; j++){
        x_prev[0][j]=0;
        x_prev[1][j]=0;
    }

    for (int i = 0; i < N; i++){
        
        for (int j = 0; j < M; j++){
	             chi = distribution(generator);
    	    x_next[1][j] = (x_prev[1][j] + dt)*r/12;
    	    x_next[0][j] = step(x_prev[0][j],x_prev[1][j],chi);
    	    x_prev[0][j] = x_next[0][j];
    	    x_prev[1][j] = x_next[1][j];
        }
     
        counter++;
        if (counter==N/1000) {
           print2file(x_prev,fptr2);
           counter = 0;
        }
     }
//    double avg = getAvgVel(x_next);
    cout << r << endl;
//    fprintf(fptr,"%15.6f\n",avg);

}
fclose(fptr2);
    return 0;
}



