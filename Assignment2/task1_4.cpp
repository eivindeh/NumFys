#include <fstream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "particle.h"
#include "int_schemes.h"
using namespace std;


int main()

{
	particle p1;
	particle p2;
	particle p3;
	const int N = 100;
	FILE* fptr;
	
	int_schemes my_schemes;
	for(int i = 0; i<N; i++){
		my_schemes::RK4(p1);
		my_schemes::Euler(p2);
		my_schemes::Mid(p3);
		
		p1::print_state(fptr);
		p2::print_state(fptr);
		p3::print_state(fptr);
	}
	
	
    return 0;
}



