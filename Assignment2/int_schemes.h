#include <armadillo>

using namespace std;
using namespace arma;

class particle;

class int_schemes {
	private:
	double dt;
	public:
	int_schemes();
	void RK4(particle *p);
	void Euler(particle *p);
	void Mid(particle *p);
	vec  midPoint(vec Bint,vec B,double dtheta);
};
