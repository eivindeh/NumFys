class particle;

class int_schemes {
	private:
	double dt;
	public:
	int_schemes();
	void RK4(particle *p);
	void Euler(particle *p);
	void Mid(particle *p);
};
