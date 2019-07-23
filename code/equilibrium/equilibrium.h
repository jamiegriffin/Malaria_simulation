
class ODE;
struct Parms;
struct Malaria_state;

struct Func_hill : public unary_function<double, double>{
	double x0;
	double x1;
	double I0;
	double k;
	Func_hill(const double x00, const double x11, const double I00, const double k0) : x0(x00), x1(x11), I0(I00), k(k0) { }
	double operator()(const double I){
		return (I==0 ? x0 : x0*((1-x1)/(1+pow(I/I0, k)) + x1));
	}
};


