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
