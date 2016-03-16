#include <iostream>
#include "int_schemes.h"
#include "particle.h"
#include <armadillo>

using namespace std;
using namespace arma;

int_schemes::int_schemes(){
	dt = 0.01;
}
void int_schemes::Euler(particle *p){
	vec x_p    = p->get_previous_state();
	vec F      = p->RHS();
	vec value  = x_p + F*dt;
	p->set_current_state(value);
	p->set_time(dt);
	
}
void int_schemes::RK4(particle *p) {
	vec k1, k2, k3, k4, x_p, x_n, F, value;
	
	x_n    = p->get_previous_state();
	F      = p->RHS();
	k1     = dt*F;
	p->set_current_state(x_n+k1/2);
	
	p->set_time(dt/2);
	
	x_p    = p->get_previous_state();
	F      = p->RHS();
	k2     = dt*F;
	p->set_current_state(x_p+k2/2);
	
	x_p    = p->get_previous_state();
	F      = p->RHS();
	k3     = dt*F;
	p->set_current_state(x_p+k3);
	
	p->set_time(dt/2);
	
	x_p    = p->get_previous_state();
	F      = p->RHS();
	k4     = dt*F;
	
	value = x_n + k1/6 + k2/3 + k3/3 + k4/6;
	p->set_current_state(value);
	
}

void int_schemes::Mid(particle *p){

	vec x_n, F, value;	
	
	x_n    = p->get_previous_state();
	F      = p->RHS();
	
	p->set_time(dt/2);
	p->set_current_state(x_n+dt/2*F);
	
	value 	= x_n + dt*p->RHS();
	
	p->set_current_state(value);
	p->set_time(dt/2);
}





