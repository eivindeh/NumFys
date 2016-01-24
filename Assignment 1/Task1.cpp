#include <fstream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iostream>
#include "functions.h"

using namespace std;

double k_bT     =  26;
double delta_U  =  80;
double alpha    =  0.2;
double L        =  20;
double tau      =  0;
double dt       =  0.01;
double r        =  12;
double eta      =  1;
int N           =  tau/dt;

double gamma_i  =  6*3.1415*r*eta;
double D_hat    =  k_bT*delta_U;
double omega    =  delta_U/(gamma_i*L*L);

void U_r(double x_hat,double t_hat,double * U)

{   double x_mod=x_hat-floor(x_hat/L)*L;
    if (x_mod > 0 && x_mod < alpha)
        *U = x_mod*k_bT/(alpha*D_hat);
    else
        *U = (1-x_mod)/(1-alpha)*k_bT/D_hat;
    
    if (t_hat > 0 && t_hat < 3*tau*omega/4)
        *U = 0;
      
}

void F_r(double x_hat,double t_hat,double * F)

{   
    double x_mod=x_hat-floor(x_hat/L)*L;

    if (x_mod > 0 && x_mod < alpha)
        *F = -k_bT/(alpha*D_hat);
    else
        *F = 1/(1-alpha)*k_bT/D_hat;
    
    if (t_hat > 0 && t_hat < 3*tau*omega/4)
        *F = 0;
      
}

double getRandom(double mu, double sigma)
{
	const double epsilon = std::numeric_limits<double>::min();
	const double two_pi = 2.0*3.14159265358979323846;

	static double z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)
	   return z1 * sigma + mu;

	double u1, u2;
	do
	 {
	   u1 = rand() * (1.0 / RAND_MAX);
	   u2 = rand() * (1.0 / RAND_MAX);
	 }
	while ( u1 <= epsilon );

	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}

double step(double x_n,double t_n)
{
    double F;
    F_r(x_n,t_n,&F);

    double chi = getRandom(1,1);

    double x_next = x_n + F*dt + sqrt(2*D_hat*dt)*chi;
    
    return x_next;
}






int main()

{
    double x[3][N];
    x[0][0]=0;
    x[1][0]=0;
    x[2][0]=0;
    
    for (int i = 0; i < N; i++){
    	x[1][i+1] = x[1][i] + dt;
    	x[0][i+1] = step(x[0][i],x[1][i]);
    	U_r(x[0][i],x[1][i],&x[2][i+1]);
    }
    cout << "hello world\n";

    return 0;
}


