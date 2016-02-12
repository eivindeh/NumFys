#include "int_schemes.h"
#include "particle.h"

int_schemes::int_schemes(){
	dt = 0.01;
}
int_schemes::Euler(particle p){
	double x_p    = p::get_previous_state()
	double F      = p::RHS()
	double value  = x_p + F*dt;
	p::set_current_state(value);
	p::update_time(dt);
	
}
