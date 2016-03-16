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
        vec v0 = {1,0,1};
        
	particle p1(1,v0,3);

	const int N = 10000;
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
