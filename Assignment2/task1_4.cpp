#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "particle.h"
#include "int_schemes.h"
#include <armadillo>

using namespace std;
using namespace arma;

int main()

{
	vec x0 = {0,0.4,0};
        vec v0 = {0.02,0,0.02};
        
	particle p1(1,x0,v0,4);

	const int N = 100000;
	FILE* fptr;
	fptr = fopen("EulerOut.txt","w");

	int_schemes my_schemes;
	for(int i = 0; i<N; i++){
    		my_schemes.RK4(&p1);
		p1.print_state(fptr);
	}
	fclose(fptr);
	
    return 0;
}
