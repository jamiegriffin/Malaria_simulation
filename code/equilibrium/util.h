//***************************************************************
//***************************************************************

struct Base_state{
	typedef double value_type;
	typedef double* iterator;
	typedef const double* const_iterator;
	typedef double& reference;
	typedef const double& const_reference;

	iterator begin(){return (iterator)this;}
	const_iterator begin() const{return (const_iterator)this;}	
	double& operator[](int i){return *(begin()+i);}
	const double& operator[](int i) const{return *(begin()+i);}
};

//***************************************************************
//***************************************************************

template<class StateT>
void set_to_zero(StateT &y){
	for(int i=0; i<y.size(); i++)
		y[i]=0;
}
//***************************************************************
//***************************************************************
template<class StateT>
void output_state(const StateT &y, ostream &s){
	const double * p =reinterpret_cast<const double *>(&y);
	for(int i=0; i<y.size()-1; i++)
		s << *p++ << '\t';
	s << *p++;
}

//***************************************************************
//***************************************************************
template<class T>
ostream& operator << (ostream &s, const pair<T, T> & x){
	s << x.first << '\t' << x.second;
	return s;
}

//***************************************************************
//***************************************************************

class Lag_state_fixed{
protected:
	double alpha;
	vector<double> history;
	int i;
	bool interpolate;
public:
	Lag_state_fixed(const double lag_length, const double update_step, const double init=0) : alpha(0), interpolate(false){
		setup(lag_length, update_step, init);
	}
	Lag_state_fixed() : alpha(0), interpolate(false) { }
	void save_state(ofstream &out){
		output_binary(i, out);
		::save_state(history, out);
	}
	void use_state(ifstream &in){
		input_binary(i, in);
		::use_state(history, in);
	}
	void clear(){
		history.clear();
		i=0;
	}
	void setup(const double lag_length, const double update_step, const double init=0){
		const int n=ceil(lag_length/update_step);
		interpolate=n*update_step-lag_length > 0.01*update_step;
		if(interpolate)
			alpha=min((n*update_step-lag_length)/update_step, 1.0);
		history.assign(n, init);
		i=0;
	}
	inline double operator()(const double x1){
		// put current value into vector and return lagged value
		const int n=history.size();
		history[i]=x1;
		i=(i+1) % n;
		return (interpolate ?  (1-alpha)*history[i] + alpha*history[(i+1) % n] : history[i]);
	}
};

//***************************************************************
//***************************************************************

class Lag_state{
protected:
public:
	double lag_length;
	double init_value;
	double init_time;
	deque<pair<double, double> > history;

	double interpolate(const double lagged_time, const int i) const{
		return 
			history[i-1].second+(lagged_time-history[i-1].first)*(history[i].second-history[i-1].second)/(history[i].first-history[i-1].first);
	}
	double lagged_value(const double lagged_time) const{
		if(history.size() < 2 || history.front().first > lagged_time){
			if(history.size()==0)
				return init_value;
			else
				return history.front().second;
		}
		int i=1;
		while(history[i].first <= lagged_time){
			i++;
			if(i == history.size())
				return history.back().second;
		}
		return interpolate(lagged_time, i);
	}
public:
	Lag_state(const double lag, const double init=0, const double init_t=0) : lag_length(lag), init_value(init), init_time(init_t) { }
	Lag_state() : lag_length(0), init_value(0), init_time(0) { }

	void save_state(ofstream &out){
		const int n=history.size();
		output_binary(n, out);
		for(int i=0; i<n; i++){
			output_binary(history[i].first, out);
			output_binary(history[i].second, out);
		}
	}
	void use_state(ifstream &in){
		int n;
		input_binary(n, in);
		history.resize(n);
		for(int i=0; i<n; i++){
			input_binary(history[i].first, in);
			input_binary(history[i].second, in);
		}
	}
	bool covers_lag(const double t) const{
		return history.size()>0 && history.front().first < t-lag_length;
	}
	void clear(){
		history.clear();
	}
	void rotate(const double t){
		for(int i=0; i<history.size(); i++)
			history[i].first += t;
	}
	void setlag(const double lag){
		lag_length=lag;
	}
	void setup(const double lag, const double init=0, const double init_t=0){
		init_value=init;
		init_time=init_t;
		history.clear();
		setlag(lag);
	}
	void update(const double current_time, const double current_value){
		if(history.size()>0 && abs(current_time-history.back().first)<1E-16)
			history.pop_back();
		history.push_back(pair<double, double>(current_time, current_value));
		while(history.size()>2 && history[1].first < current_time-lag_length)
			history.pop_front();
	}
	double operator()(const double current_time) const{
		return lagged_value(current_time-lag_length);
	}
	double operator()(const double current_time, const double lag_length_b) const{
		if(lag_length_b > lag_length)
			error_crit("lag_length_b must not be greater than lag_length in Lag_state::()(current_time, lag_length_b)");
		return lagged_value(current_time-lag_length_b);
	}
	double mean(const double u0){
		if(history.size()==0 || history.front().first > u0)
			error_crit("History doesn't cover time span in Lag_state::mean");

		const double eps=1E-6;
		int i0=1;
		while(history[i0].first <= u0-eps)
			i0++;
		
		double S=0;
		double y0=interpolate(u0, i0);
		double t0=u0;
		for(int i=i0; i<history.size(); i++){
			double t1=history[i].first;
			double y1=history[i].second;	
			S += 0.5*(y1+y0)*(t1-t0);
			t0=t1;
			y0=y1;
		}
		return S/(t0-u0);
	}
};

//***************************************************************
//***************************************************************

template<class T>
T max(const T &a1, const T &a2, const T &a3){
	return max(max(a1, a2), a3);
}
template<class T>
T max(const T &a1, const T &a2, const T &a3, const T &a4){
	return max(max(a1, a2), max(a3, a4));
}
template<class T>
T max(const T &a1, const T &a2, const T &a3, const T &a4, const T &a5){
	return max(max(a1, a2), max(a3, a4, a5));
}
template<class T>
T max(const T &a1, const T &a2, const T &a3, const T &a4, const T &a5, const T &a6){
	return max(max(a1, a2, a3), max(a4, a5, a6));
}
template<class T>
T max(const T &a1, const T &a2, const T &a3, const T &a4, const T &a5, const T &a6, const T &a7){
	return max(max(a1, a2, a3), max(a4, a5, a6, a7));
}
template<class T>
T max(const T &a1, const T &a2, const T &a3, const T &a4, const T &a5, const T &a6, const T &a7, const T &a8){
	return max(max(a1, a2, a3, a4), max(a5, a6, a7, a8));
}

template<class T>
T max_abs(const T &a1, const T &a2){
	return max(abs(a1), abs(a2));
}
template<class T>
T max_abs(const T &a1, const T &a2, const T &a3){
	return max_abs(max_abs(a1, a2), a3);
}
template<class T>
T max_abs(const T &a1, const T &a2, const T &a3, const T &a4){
	return max_abs(max_abs(a1, a2), max_abs(a3, a4));
}
template<class T>
T max_abs(const T &a1, const T &a2, const T &a3, const T &a4, const T &a5){
	return max_abs(max_abs(a1, a2), max_abs(a3, a4, a5));
}
template<class T>
T max_abs(const T &a1, const T &a2, const T &a3, const T &a4, const T &a5, const T &a6){
	return max_abs(max_abs(a1, a2), max_abs(a3, a4, a5, a6));
}
template<class T>
T max_abs(const T &a1, const T &a2, const T &a3, const T &a4, const T &a5, const T &a6, const T &a7){
	return max_abs(max_abs(a1, a2, a3), max_abs(a4, a5, a6, a7));
}
template<class T>
T max_abs(const T &a1, const T &a2, const T &a3, const T &a4, const T &a5, const T &a6, const T &a7, const T &a8){
	return max_abs(max_abs(a1, a2, a3), max_abs(a4, a5, a6, a7, a8));
}

//***************************************************************
//***************************************************************

struct Lookup_table{
	double min_x;
	double max_x;
	double dx;
	
	double* y;
	double* z;
	double r;
	int n;
	bool interpolate;

	Lookup_table() : interpolate(true), y(0), n(0) {}
	~Lookup_table(){
		reallocate(0);
	}

	Lookup_table& operator=(const Lookup_table &t){
		reallocate(t.n);
		min_x=t.min_x;
		max_x=t.max_x;
		dx=t.dx;
		r=t.r;
		interpolate=t.interpolate;
		copy(t.y, t.y+n, y);
		if(n>1)
			copy(t.z, t.z+n-1, z);
		return *this;
	}
	void reallocate(const int n1=0){
		if(n>0)
			delete [] y;
		if(n>1)
			delete [] z;
		n=n1;
		if(n>0)
			y=new double [n];
		if(n>1)
			z=new double [n-1];
	}
	void setup(const string &file){
		interpolate=true;
		vector<double> x, y0;
		import_file(file, x, y0);
		n=y0.size();
		reallocate(n);
		
		min_x=x.front();
		max_x=x.back();
		dx=(max_x-min_x)/(x.size()-1);
		r=1/dx;
		if(dx<=0)
			error_crit("First input column should be increasing in Lookup_table::setup, input from file "+file);
		for(int i=1; i<x.size(); i++)
			if(abs(x[i]-x[i-1] - dx)/dx>1E-2)
				error_crit("Step size in first input column should be constant in Lookup_table::setup, input from file "+file);
		copy(y0.begin(), y0.end(), y);
		for(int i=0; i<n-1; i++)
			z[i]=(y[i+1]-y[i])*r;
	}
	template<class F>
	void setup(F f, const double n0, const double min_x0, const double max_x0){
		interpolate=true;
		reallocate(n0+1);
		min_x=min_x0;
		max_x=max_x0;
		dx=(max_x-min_x)/n0;
		r=1/dx;
		for(int i=0; i<=n0; i++)
			y[i]=f(min_x + i*dx);
		for(int i=0; i<n0; i++)
			z[i]=(y[i+1]-y[i])*r;
	}

	inline double operator()(const double x) const{
		if(x<=min_x)
			return y[0];
		if(x>=max_x)
			return y[n-1];

		int i=(x-min_x)*r;		
		if(interpolate){
			const int n1=n-2;
			i=i>n1 ? n1 : i;
			const double x0=min_x+dx*i;
			return y[i] + (x-x0)*z[i];
		}
		else{
			const int n1=n-1;
			i=i>n1 ? n1 : i;
			return y[i];
		}
	}
};

struct Pow_func : public unary_function<double, double>{
	double a;
	Pow_func(const double a0) : a(a0) { }
	double operator()(const double x){
		return pow(x, a);
	}
};
//***************************************************************
//***************************************************************

