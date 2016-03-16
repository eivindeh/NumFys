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
	current_state = zeros<vec>(6);
	previous_state =zeros<vec>(6);
	t=0;
}

particle::particle(double _q, vec v0, int _task){ 
	q = _q;
	task = _task;
	current_state = zeros<vec>(6);
	previous_state = zeros<vec>(6);
        previous_state(span(3,5))= v0;
        setField();
	t=0;
}

void particle::setField(){
	double B0;
	double beta,x,y;
	mat R;
	vec B_cyl;
	switch(task){
		case 1:
			B  = {0,0,1};
			E  = {0,0.1,0.1};
		break;
		
		case 2:	
			B0   = 1;
			beta = 0.1;
		
			B    = {0,0,B0 + beta*previous_state(4)};
			E    = {0,0,0};
		break;
		
		case 3:
			x = previous_state(2);
			y = previous_state(3);
			B0 = 1;
			R     = {{x/sqrt(x*x+y*y), -y/sqrt(x*x+y*y), 0},
				 {y/sqrt(x*x+y*y), x/sqrt(x*x+y*y),  0},
				 {0		 , 0               , 1}};
				 
			B_cyl = {0, B0, 0};
			B = R*B_cyl;
		
			E = {0, 0, 0};
		break;
		
		default:
		cout<<"something is wrong"<<endl;
		
	}	
	
}

vec particle::RHS(){
        setField();
        vec F = zeros<vec>(6);
        F(span(0,2))=previous_state(span(3,5));
        F(span(3,5))=(E+cross(previous_state(span(3,5)),B));
	/*mat A = zeros<mat>(6,6);
        A(span(0,2),span(3,5)) = eye<mat>(3,3);
	A(3,4)=q*B;
	A(4,3)=-q*B;*/
	return F;
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


