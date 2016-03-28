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
//	vec x0 = {0,0.4,0};
//      vec v0 = {0.022,0,0.02};

	vec x0 = {0,4,0};
        vec v0 = {1,0,1};
        v0 = v0*0.1/sqrt(2);
        
        int task = 5;
        // Task number specifies which task to perform
	// Task 1 = constant B ans E
	// Task 2 = zero E and B(y)
	// Task 3 = zero E and B(r)
	// Task 4 = Magnetic bottle
	// Task 5 = Earth field
	// Task 6 = Earth field whith guiding center approximation
	  
        
	particle p1(-1,x0,v0,task);
	const int N = 3000000;
	FILE* fptr;
	fptr = fopen("EulerOut.txt","w");

	int_schemes my_schemes;
	for(int i = 0; i<N; i++){
    		my_schemes.RK4(&p1);
		p1.print_state(fptr,100);
		
		if(i%10000==0){  
        	cout<<i*100.0/N<<" %"<<endl;
        	}
	}
	fclose(fptr);
	
    return 0;
}
