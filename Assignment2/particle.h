#include <armadillo>
using namespace std;
using namespace arma;

class particle {
	private:
	vec 	B;
	vec	E;
	double 	m;
	int     n;
	double  mu;
	vec     v;
	double 	q;
	double 	t;
	double  j;
	vec 	current_state;
	vec 	previous_state;
	public:
	particle();
	particle(double _q,vec x0, vec v0, int task);
	int     task;
	vec 	RHS();
	vec     return_integrand(double theta);
	void 	print_state(FILE* fptr,int M);
	void 	set_current_state(vec new_current_state);
	vec 	get_previous_state();
	void 	set_time(double time);
	void	setField();
}
;
