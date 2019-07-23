
//***************************************************************
//***************************************************************

struct Mosq_state{
	double& S;
	double& E;
	double& I;
	double& EL;
	double& LL;
	double& PL;
	Mosq_state(double *y) : S(*y), E(*(y+1)), I(*(y+2)), EL(*(y+3)), LL(*(y+4)), PL(*(y+5))  {
		set_init();
	}
	void operator = (const Mosq_state &z){}
	void set_init(){
		S=1;
		E=0;
		I=0;
		EL=1;
		LL=1;
		PL=1;
	}
	static int size(){
		return 6;
	}
};

struct Human_state{
	double* S;	
	double* T;	
	double* D;	
	double* A1;	
	double* U;
	double* P;
	double* IC_A;
	double* IB;	
	double* IV;	
	double* ID;	

	int n_age;

	int n_prev_states;
	int n_states;

	Human_state(const int na=-1) : n_age(na) {
		n_prev_states= 6;
		n_states=n_prev_states + 4;
	}
	void setup(double *y){
		for(int l=0; l<size(); l++)
			y[l]=0;

		int k=0;
		S=	y + (k++)*n_age;
		T=	y + (k++)*n_age;
		D=	y + (k++)*n_age;
		A1=	y + (k++)*n_age;
		U=	y + (k++)*n_age;
		P=	y + (k++)*n_age;
		assert(k==n_prev_states);
		IC_A=y+ (k++)*n_age;
		IB=	y + (k++)*n_age;
		IV=	y + (k++)*n_age;
		ID=	y + (k++)*n_age;
		assert(k==n_states);
		assert(k*n_age==size());
	} 
	inline int size() const{
		return n_states*n_age;
	}
	inline int num_states() const{
		return n_states;
	}
	inline int num_prev_states() const{
		return n_prev_states;
	}
	
	inline double* state(const int k){
		return S + k*n_age;
	}
	inline const double* state(const int k) const{
		return S + k*n_age;
	}

	double& operator[](const int l){
		return S[l];
	}
	const double& operator[](const int l) const{
		return S[l];
	}
};

struct Malaria_state{

	double* y;

	vector<Mosq_state> mosq_state;
	vector<vector<Human_state> > human_state;
	int num_het;
	int num_age;
	int total_size;
	int current_size;
	int init_size;
	Malaria_state() : num_het(-1), num_age(-1), y(0) {	}
	~Malaria_state() {
		delete_storage();
	}
	void delete_storage();
	void setup(const Malaria_state &z);
	void setup(const int na, const int nh);
	void setup(const Parms &parms);
	
	void add_intervention(const double coverage);
	
	void operator = (const Malaria_state &z);

	inline int size() const{		
		return total_size; 
	}
	inline int size_now() const{		
		return current_size; 
	}
	inline double& operator[](const int l){
		return y[l];
	}
	inline const double& operator[](const int l) const{
		return y[l];
	}
};

//***************************************************************
//***************************************************************