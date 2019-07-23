//******************************************************************************
//******************************************************************************

struct Mosq_pop{

	void save_state(ofstream &out);
	void use_state(ifstream &in);

	//***************************************************************************
	// needed for deterministic solving:
	Parms* parms;
	double infv;
	Lag_state lag_infv;
	
	Ran ran;

	void update_history(const double t);

	int N;

	int init_pop_year;
	vector<double> pop_r;
	bool growing_pop;

	vector<double> b_rate;
	vector<double> b_rate_t;
	int birth_index;
	
	//***************************************************************************

	void change_itn(const int n);
	void change_irs(const int n);

	Village* village;

	string species;

	// fixed parameters
	double latgam;
	double latmosq;
	double omega;
	bool stochastic;
	bool feedback;
	double itn_decay_rate;
	double irs_decay_rate;
	double itn_decay_shape;
	double irs_decay_shape;

	double mue;
	double mul;
	double de;
	double dl;
	double dp;
	double eov;
	double mup;
	double gamma;
	double b_lambda;

	// natural life cycle parameters
	double mu0;
	double tau1;
	double tau2;

	vector<double> seasonal_table;
	double seasonal_curve(const double t) const;

	void find_beta();
	void find_init();

	double find_eir(const double t);

	//***************************************************************************
	double beta;
	
	bool first;
	double FOIM;
	double FOIM_end;
	Lag_state_fixed lag_FOIM;
	Lag_state_fixed lag_eir;		
	Lag_state_fixed lag_inc_E;

	
	double FOIM_init;
	double C;
	double rel_M;

	double K0;
	double M0;
	
	double S;
	double E;
	double I;

	double EL;
	double LL;
	double PL;

	double p1;
	double p2;

	//***************************************************************************
	// larval habitat reduction
	bool larvicide;
	double larval_time;
	double deprec_param;
	double larvi_min;
	void change_larvicide(const double p, const int num);

	//***************************************************************************

	// life cycle parameters possibly changed by ITN, IRS
	double av;
	double mu;
	double PM;

	double f;
	double Q;

	double D;
	double w_bar;
	double z_bar;
	//***************************************************************************
	double Q0;
	double Q_in;
	double Q_bed;

	double itn_decay(const double t_ago);
	double irs_decay(const double t_ago);

	
	bool itn_irs_2016;
	bool irs_ellie;

	int year0;

	double itn_rn;
	double itn_rnw;
	double itn_dnw;
	double itn_dnf;
	double irs_ri;
	double irs_diw;
	double irs_dif;
	double irs_riw;

	double endophily;
	
	double itn_kill;
	double itn_repel;
	double itn_repel_min;
	double irs_kill;
	double irs_repel;

	double irs_decay_mort1;
	double irs_decay_mort2;
	double irs_decay_succ1;
	double irs_decay_succ2;
	double irs_decay_det1;
	double irs_decay_det2;
	double irs_k0;
	bool irs_succ_k0;

	double prob_fail(const ITN_IRS &itn_irs, const double t);
	double prob_repeat(const ITN_IRS &itn_irs, const double t);
	double prob_bite(const ITN_IRS &itn_irs, const double t);

	Mosq_pop(Village* v=0);
	void update(const double dt);
	void update_foim(const double sum_fail, const double sum_repeat, const double sum_rel_biting_rate);
	void setup(Parms* parms1, const string &name, const int N=1, const string &pop_file="");
	void reset();
	void set_init_indiv(Mosq_state &mosq_state, const double t);
	
	void update_det(const Mosq_state &y,  Mosq_state &dydt, const double t);
	void find_equilibrium_det(Mosq_state &y_mosq, const double mean_c);
};

//******************************************************************************
//******************************************************************************