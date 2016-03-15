#include <armadillo>
using namespace std;
using namespace arma;

class particle {
	private:
	double 	B;
	double 	m;
	double 	v0_x;
	double 	v0_y;
	double 	v0_z;
	double 	q;
	double 	t;
	vec 	current_state;
	vec 	previous_state;
	public:
	particle();
	particle(double _B, double _q, double v0_x,double m, double v0_z);

	vec 	RHS();
	void 	print_state(FILE* fptr);
	void 	set_current_state(vec new_current_state);
	vec 	get_previous_state();
	void 	set_time(double time);
}
;
