
class particle {
	private:
	double B;
	double v0_x;
	double v0_y;
	double v0_z;
	double q;
	double current_state[6];
	double previos_state[6];
	double t;
	public:
	particle();
	particle(double _B, double _q, double v0_x, double v0_z);
	double* RHS();
	void print_state(FILE* fptr);
	void set_current_state(double new_current_state[]);
	double* get_previous_state();
	void set_time(double time);
};
