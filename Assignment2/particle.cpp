#include "particle.h"
#include <array>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>

particle::particle(){
    B=1;
	q=1;
	m=1;
	current_state={0};
	previous_state={0};
	t=0;
}

particle::particle(double _B, double _q, double v0_x, double v0_z){ 
	B = _B;
	q = _q;
	current_state={0};
	previous_state={0,0,0,v0_x,0,v0_z};
	t=0;
}

double* particle::RHS(){
	double f[6];
	f = {previous_state[4],
		 previous_state[5],
		 previous_state[6],
		 previous_state[5]*q*B/m,
		-previous_state[4]*q*B/m,
		 0};
	return f;
}

void particle::print_state(FILE *fptr){  
        fprintf(fptr,"%15.6f%15.6f%15.6f\n",current_state[1],current_state[2],current_state[3]);
}

void particle::set_current_state(double new_current_state[]){
		previous_state=new_current_state ;
	}

double* particle::get_previous_state(){
	return previous_state;
}

void update_time(double dt){
	t = t+dt;
}
