#include "particle.h"
#include "int_schemes.h"
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
	//variables used in case 2
	double B0, beta;
	
	//variables used in case 3
	double x,y;
	mat R;
	vec B_cyl;
	
	//variables used in case 4
	double theta_max = 2*3.1415;
	int M = 100;
	double dtheta = theta_max/M;
	double theta = 0;
	int i;
	int_schemes my_scheme;
	
	switch(task){
		//Set constant B and E field
		case 1:
			B  = {0,0,1};
			E  = {0,0.1,0.1};
			break;
		
		//Set zero E field and B = B(y) e_y		
		case 2:	
			B0   = 1;
			beta = 0.1;
		
			B    = {0,0,B0 + beta*previous_state(4)};
			E    = {0,0,0};
			break;
		
		//Set zeros E field and B = B(r) e_theta
		case 3:
			x = previous_state(0)+0.000001;
			y = previous_state(1)+0.000001;
			B0 = 0.5;
			R     = {{x/sqrt(x*x+y*y), -y/sqrt(x*x+y*y), 0},
				 {y/sqrt(x*x+y*y), x/sqrt(x*x+y*y),  0},
				 {0		 , 0               , 1}};
				 
			B_cyl = {0, B0, 0};
			B = R*B_cyl;
		
			E = {0, 0, 0};
			break;
		
		//Set helmholtz field 
		case 4:
			B = zeros<vec>(3);
			x = previous_state(0)+0.000001;
			y = previous_state(1)+0.000001;
			for (i = 0; i<M; i++){
				B = my_scheme.midPoint(return_integrand(theta),B,dtheta);
				theta = theta+dtheta;
			}
			
			R     = {{x/sqrt(x*x+y*y), -y/sqrt(x*x+y*y), 0},
				 {y/sqrt(x*x+y*y), x/sqrt(x*x+y*y),  0},
				 {0		 , 0               , 1}};
				 
		        B = R*B;
			
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

vec particle::return_integrand(double theta){
	//only in use in case 4
	double r 	= norm(previous_state(span(0,1)),2);
	double R	= 1;
	double z 	= previous_state(2);
	
	vec B_int = zeros<vec>(3);
	
	B_int(0) = pow(1/(R*R),3/2)/(4*3.1415*R)*((z-1)*cos(theta)/pow((pow((r-R*cos(theta)),2)+pow(R*sin(theta),2)+pow((z-1),2)),(3/2)) + (z+1)*cos(theta)/pow((pow((r-R*cos(theta)),2)+pow(R*sin(theta),2)+pow((z+1),2)),(3/2))); 
	
	B_int(2) = pow(1/(R*R),3/2)/(4*3.1415*R)*((1-cos(theta)*r/R)/pow((pow((r-R*cos(theta)),2)+pow(R*sin(theta),2)+pow((z-1),2)),(3/2)) + (1-cos(theta)*r/R)/pow((pow((r-R*cos(theta)),2)+pow(R*sin(theta),2)+pow((z+1),2)),(3/2)));
	
	return B_int;
}

