#include "particle.h"
#include <array>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

particle::particle(){
        B=1;
	q=1;
	m=1;
	current_state = zeros<vec>(6);
	previous_state =zeros<vec>(6);
	t=0;
}

particle::particle(double _B, double _q, double _m, double v0_x, double v0_z){ 
	B = _B;
	q = _q;
	m = _m;
	current_state = zeros<vec>(6);
	previous_state = zeros<vec>(6);
        previous_state(3) = v0_x;
	previous_state(5) = v0_z;
	t=0;
}

vec particle::RHS(){
	mat A = zeros<mat>(6,6);
        A(span(0,2),span(3,5)) = eye<mat>(3,3);
	A(3,4)=q*B/m;
	A(4,3)=-q*B/m;
	return A*previous_state;
}

void particle::print_state(FILE *fptr){  
        fprintf(fptr,"%15.6f%15.6f%15.6f\n",previous_state(0),previous_state(1),previous_state(2));
}

void particle::set_current_state(vec new_current_state){
	previous_state=new_current_state;
	}

vec particle::get_previous_state(){
	return previous_state;
}

void particle::set_time(double dt){
	t = t+dt;
}


