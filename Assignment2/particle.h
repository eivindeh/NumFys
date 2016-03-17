#include <armadillo>
using namespace std;
using namespace arma;

class particle {
	private:
	vec 	B;
	vec	E;
	double 	m;
	vec     v;
	double 	q;
	double 	t;
	int     task;
	vec 	current_state;
	vec 	previous_state;
	public:
	particle();
	particle(double _q,vec x0, vec v0, int task);

	vec 	RHS();
	vec     return_integrand(double theta);
	void 	print_state(FILE* fptr);
	void 	set_current_state(vec new_current_state);
	vec 	get_previous_state();
	void 	set_time(double time);
	void	setField();
}
;
