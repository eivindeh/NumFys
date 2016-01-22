#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void U_r();

int main(argc, argv)
int argc;
char *argv[];
{
    double U =0;
    U_r(&U);
    printf("%f",U);

    return 0;
}

void U_r(double * U)

{
    *U=3;
}
