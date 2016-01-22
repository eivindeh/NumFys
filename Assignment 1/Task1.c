#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double k_b=2;
double T=6;
double alpha=0.4;
double L=2;
double D_hat=3;
double tau=10;
double omega=2;

void U_r();
void F_r();

int main(argc, argv)
int argc;
char *argv[];
{
    double U =0;
    double x_hat=3;
    double t_hat=40;
    F_r(x_hat,t_hat,&U);
    printf("Hello world %f \n",U);

    return 0;
}

void U_r(double x_hat,double t_hat,double * U)

{   
    if (x_hat > 0 && x_hat < alpha)
        *U = x_hat*k_b*T/(alpha*D_hat);
    else
        *U = (1-x_hat)/(1-alpha)*k_b*T/D_hat;
    
    if (t_hat > 0 && t_hat < 3*tau*omega/4)
        *U = 0;
      
}

void F_r(double x_hat,double t_hat,double * F)

{   
    if (x_hat > 0 && x_hat < alpha)
        *F = -k_b*T/(alpha*D_hat);
    else
        *F = 1/(1-alpha)*k_b*T/D_hat;
    
    if (t_hat > 0 && t_hat < 3*tau*omega/4)
        *F = 0;
      
}
