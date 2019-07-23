template<class F>
struct F_map{
	F &f;
	int n_tries;
	map<double, double> previous_f;
	bool bracketed;
	F_map(F &f0) : f(f0), n_tries(0), bracketed(false) { }  
	double operator()(const double x){
		if(n_tries++==100)
			error_crit(string("Too many function calls in find_root") 
					   + (bracketed ? string(" after bracketing root") : string(" while trying to bracket root")));
		map<double, double>::const_iterator it=previous_f.find(x);
		return it!=previous_f.end() ? it->second : previous_f[x] = f(x);
	}
};

template<class F>
double find_root(F &f, const double x_guess, const double tol, const double r0, const double r1){
	// find the root of an increasing function f(x), given an initial guess x_guess
	
	double x1=x_guess;
	double x0=x1;
	F_map<F> f_map(f);
	
	if(abs(f_map(x1))<tol)
		return x1;

	while(f_map(x0) > 0){
		x1=x0;
		x0 -= r0;
	}
	while(f_map(x1) < 0){
		x0=x1;
		x1 += r1;
	}
	f_map.bracketed=true;
	
	double x2, f2;
	do{
		const double f0 = f_map(x0);
		const double f1 = f_map(x1);
		x2 = (x0*f1 - x1*f0)/(f1-f0);
		f2 = f_map(x2);
		if(f2>0) x1=x2;					
		else x0=x2;
	}while(abs(f2)>tol);
	return x2;
}
