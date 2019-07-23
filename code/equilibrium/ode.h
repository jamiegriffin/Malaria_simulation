//***************************************************************
//***************************************************************

class ODE{
public:
	void setup(Malaria_state &y0, Parms& parms0, const double t0, const double h0);
	void step(const double total_step);
	double time() const {
		return t;
	}
public:
	int step_type;
	double t;
	double step_size;
	Malaria_state ytemp;
	Malaria_state k1;
	Malaria_state k2;
	Malaria_state k3;
	Malaria_state k4;

	Malaria_state* y;
	Parms* parms;

	void resize();
	void sub_step(const double h1);
};

//***************************************************************
//***************************************************************
