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

mat getTransMat(particle p){
	double x,y,z,sf,cf,st,ct,r;
	mat R;
	vec X = p.get_previous_state();
	x  = X(0);
	y  = X(1);
	z  = X(2);
	r = norm(X(span(0,2)),2);
	if(p.task == 5 || p.task == 6)
	{
		sf = y/sqrt(x*x+y*y);
		cf = x/sqrt(x*x+y*y);
		st = sqrt(r*r-z*z)/r;
		ct = z/r; 
		R = {{st*cf, ct*cf, -sf},
		     {st*sf, ct*sf,  cf},
		     {ct   ,-st   ,  0}};
	     }
	else if(p.task ==4 || p.task == 3){
		R     = {{x/sqrt(x*x+y*y), -y/sqrt(x*x+y*y), 0},
			 {y/sqrt(x*x+y*y), x/sqrt(x*x+y*y),  0},
			 {0		 , 0               , 1}};
		}
	return R;
}

particle::particle(){
        B=1;
	q=1;
	current_state = zeros<vec>(6);
	previous_state =zeros<vec>(6);
	t=0;
}

particle::particle(double _q,vec x0, vec v0, int _task){ 
	q = _q;
	j = 0;
	n = 0;
	mu = 0;
	task = _task;
	current_state = zeros<vec>(6);
	previous_state = zeros<vec>(6);
        previous_state(span(3,5))= v0;
        previous_state(span(0,2))= x0;
        setField();
	t=0;
}

void particle::setField(){
	//variables used in case 2
	double B0, beta;
	
	//variables used in case 3
	mat R;
	vec B_cyl;
	
	//variables used in case 4
	double theta_max = 2*3.1415;
	int M = 10;
	double dtheta = theta_max/M;
	double theta = 0;
	int i;
	int_schemes my_scheme;
	
	//case 5
	double r, r3;
	vec B_s;
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
			B0 = 1;
			R = getTransMat(*this);
				 
			B_cyl = {0, B0, 0};
			B = R*B_cyl;
		
			E = {0, 0, 0};
			break;
		
		//Set helmholtz field 
		case 4:
			B = zeros<vec>(3);

			theta = 0;
			for (i = 0; i<M; i++){
				B = B + my_scheme.simpsons(this, theta, dtheta);
				theta = theta+dtheta;
			}
			R = getTransMat(*this);		 
		        B = R*B;
			E = {0, 0, 0};
			break;
			
		//set earth field
		case 5:
		        B_s = zeros<vec>(3);
		        B = zeros<vec>(3);
			r = norm(previous_state(span(0,2)),2);
			r3 = pow(r,3);
			theta = acos(previous_state(2)/r);
			B_s(0) = -2*cos(theta)/r3;
			B_s(1) = -sin(theta)/r3;
			R = getTransMat(*this);
		        B = R*B_s;
		        E = {0, 0, 0};
		        break;
		default:
		cout<<"something is wrong"<<endl;
		
	}	
	
}

vec particle::RHS(){
	vec F = zeros<vec>(6);
	double v_par,v_nor,r,v2,B_abs,alpha,theta,r3;
	mat R;
	vec b,gradB;
        if(task==6){
        //using guiding center approximation
		v_par = previous_state(3);
		r = norm(previous_state(span(0,2)),2);
		theta = acos(previous_state(2)/r);
		r3 = pow(r,3);
		B = zeros<vec>(3);
        	B(0) = -2*cos(theta)/r3;
		B(1) = -sin(theta)/r3;
		B(2) = 0;
		B_abs = norm(B,2);
		R = getTransMat(*this);
		b = B/B_abs;
		gradB ={-3*sqrt((4-3*pow(sin(theta),2)))/pow(r,4),
			-3*cos(theta)/(pow(r,4)*sqrt((4-3*pow(sin(theta),2)))),
			 0};
		if (n==0) {
        		v_par = dot(R*b,previous_state(span(3,5)));
        		v_nor = norm(previous_state(span(3,5))-v_par*R*b,2);
        		cout<<r<<endl;
        		mu = v_nor*v_nor/B_abs;
        		n = 1;
        		 cout<<"gradB:"<<gradB<<endl;
        	}
		v2 = v_par*v_par+mu*B_abs;
		//alpha=1;
		alpha=62.56;
		F(span(0,2))= R*(v2/(2*alpha*B_abs*B_abs)*(1+v_par*v_par/v2)*(cross(b,gradB))+v_par*b);
		F(3) = -1.0/2.0*(dot(R*b,R*gradB))*mu;
		//cout<<F(3)<<endl;
	}
	
	else{
		setField();
		F(span(0,2))=previous_state(span(3,5));
		F(span(3,5))=(E+cross(previous_state(span(3,5)),B));
		}
        
	return F;
	
}

void particle::print_state(FILE *fptr, int K){
        if(j==K){  
        	fprintf(fptr,"%15.6f%15.6f%15.6f\n",previous_state(0),previous_state(1),previous_state(2));
        	j=1;
        	}
        else{
        	j++;
        }
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
	double R	= 0.2;
	double z 	= previous_state(2);
	
	vec B_int = zeros<vec>(3);
	
	B_int(0) = pow(1+(R*R),3.0/2.0)/(4*3.1415*R)*((z-1)*cos(theta)/pow((pow((r-R*cos(theta)),2)+pow(R*sin(theta),2)+pow((z-1),2)),(3.0/2.0)) + (z+1)*cos(theta)/pow((pow((r-R*cos(theta)),2)+pow(R*sin(theta),2)+pow((z+1),2)),(3.0/2.0))); 
	
	B_int(2) = pow(1+(R*R),3.0/2.0)/(4*3.1415)*((1-cos(theta)*r/R)/pow((pow((r-R*cos(theta)),2)+pow(R*sin(theta),2)+pow((z-1),2)),(3.0/2.0)) + (1-cos(theta)*r/R)/pow((pow((r-R*cos(theta)),2)+pow(R*sin(theta),2)+pow((z+1),2)),(3.0/2.0)));
	
	return B_int;
}


