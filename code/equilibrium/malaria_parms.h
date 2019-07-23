//******************************************************************************
//******************************************************************************

struct Parms : public Base_state{

public:
	
	Malaria_state state;
	ODE ode;
	Simulation* simulation;
	
	double prev_age0;
	double prev_age1;
	Lag_state lag_prev;

	/* run_type
	0	run full model (with larval stages if feedback == 1)
	1	run transmission model, using stored mosquito birth rate
	2	run just larval model and store adult mosquito birth rate
	*/
	
	int run_type;

	void calc_dydt(const double t, Malaria_state &y, Malaria_state &dydt);
	ITN_IRS itn_irs_det;
	double itn_coverage;
	double irs_coverage;
	double update_coverage_time;
	bool leave_itn_on;
	double leave_itn_time;
	
	void change_itn_coverage(const double p1);
	void add_itn(const double coverage);
	void add_irs(const double coverage);

	//**********************************************************************
	//**********************************************************************

	void input_full();
	void update_foim(const bool update_lag, const double t);
	void find_eir(const Malaria_state &y, const double t=-36500.0);

	//**********************************************************************
	//**********************************************************************
	

	// biting rate by age: FOI(a) = FOI0*(1-rho*exp(-a/a0))
	double rho;		// 1 - intercept for FOIh at age 0
	double a0;		// parameter for increase in FOIh with age
	
	// duration of infection states
	double dur_E;	// time in latent period
	double dur_I;	// time in untreated disease plus patent asymptomatic infection (dur_I = dur_D + dur_A0)
	double dur_T;	// time in clinical disease with treatment (i.e. time with patent parasitaemia if had clinical disease and were treated)
	double dur_D;	// time in clinical disease without treatment (i.e. time from patent parasitaemia with clinical disease until end of clinical disease)
	double dur_U;	// time in subpatent infection
	double dur_P;	// time in prophylaxis

	// infection, anti-infection immunity
	double bh;		// probability a human becomes infected upon being bitten by an infectious mosquito
	double bmin;
	double db;
	double kb;
	double IB0;
	double ub;
	double gammab;
	double fb0;
	double ab0;

	// detection immunity
	double dmin;
	double dd;
	double kd;
	double ID0;
	double ud;
	double gammad;
	double fd0;
	double ad0;

	// clinical immunity
	double phi0;
	double phi1;
	double dc;		// duration parameter for acquired clinical immunity
	double kc;
	double IC0;
	double uc;
	double gammac;
	double fc0;
	double ac0;

	double P_IC_M;	// maternal clinical immunity as proportion of mother's
	double dm;		// duration parameter for maternal clinical immunity

	// severe disease immunity
	double theta0;
	double theta1;
	double dv;
	double kv;
	double IV0;
	double uv;
	double gammav;
	double fv0;
	double av0;
	double P_IV_M;
	double dvm;
	
	double tau_v;

	// variation in exposure
	double sigma2;

	double latgam;	// lag from asexual parasites to mature gametocytes, doesn't affect equilibrium (as currently programmed, conditional on FOIh)

	// Infection parameters
	// human to mosquito
	// probability that mosquito becomes infected
	double cT;	
	double cD;	
	double cU;	

	double gamma_inf;

	double alpha_pcr;
	double beta_pcr;

	//**********************************************************************
	//**********************************************************************

	vector<Mosq_pop>* mosq_pop;

	int num_age;
	int num_het;

	// Human malaria parameters
	double ft;			// proportion with clinical disease successfully treated for a fully susceptible genotype
	// parameters derived from other parameters

	double omega;		// biting rate normalisation 


	double rec_A0;		// natural recovery rate from asymptomatic patent to sub-patent parasitaemia
	double rec_T;		// recovery rate from clinical malaria with treatment
	double rec_D;		// recovery rate from clinical malaria to asymptomatic patent parasitaemia without treatment
	double rec_U;		// recovery rate from sub-patent infection
	double rec_P;		// rate of leaving prophylaxis

	
	double dec_dc;
	double dec_db;
	double dec_dd;
	double dec_dv;

	//**********************************************************************

	int inf_immunity;
	int det_immunity;
	int clin_immunity;
	int sev_immunity;

	bool separate_sev;

	// variables calculated initially
	vector<double> fB;
	vector<double> fC;
	vector<double> fD;
	vector<double> fV;

	vector<double> fB1;
	vector<double> fC1;
	vector<double> fD1;
	vector<double> fV1;

	vector<double> IC_M_age;
	vector<double> IV_M_age;
	int age_20_l;
	int age_20_u;
	double age_20_factor;

	vector<double> rel_foi;
	vector<double> het_wt;
	vector<double> het_x;
	vector<double> x_I;
	vector<double> FOIh_age; // ratio of FOIh at this age to maximum FOIh

	Lookup_table b_table;
	Lookup_table c_table;
	Lookup_table d_table;
	Lookup_table v_table;

	//**********************************************************************
	// variables derived from model output
	double Sh;	// sums across ages in each human state
	double Th;
	double Dh;
	double A1h;
	double Ph;
	double Uh;
	double H_sum;

	vector<vector<double> > IC_20;
	vector<vector<double> > IV_20;
	vector<double> clin_inc;
	vector<double> inf_inc;
	vector<double> sev_inc;
	//Human_state H_sum_age;
	
	vector<double> EIR;	// daily EIR
	vector<double> mean_c;

	// lagged state variables
	vector<Lag_state> lag_EIR;	// lagged by human latent period
	vector<Lag_state> lag_mean_c;			// lagged by delay from asexual parasites to gametocytes

	//******************************************************************************
	//******************************************************************************
	// Human demographic parameters and variables
	// demographic variables
	vector<double> b;	// widths of age categories
	vector<double> age;	// age at start of age category
	vector<double> N_pop;	// vector of proportion of population in each age group

	// rate of exiting age groups, death rate plus growth rate
	vector<double> alpha;

	Guide_table age_guide;
			
	double eta; // birth and death rate

	double birth_rate;

	vector<double> age_rate;
	mutable vector<double> temp_pred;

	//******************************************************************************
	//******************************************************************************
	Parms();
	
	double prop_in_age_range(const double age1, const double age2) const;

	void setup(const int na, const int nh, const vector<vector<double> > &v, const double birth_rate0);
	void reset();

	void update_history(Malaria_state &y, const double current_time);
	void clear_history();
	void rotate_history(const double t);
	void copy_history(const Parms &parms);

	void update_parms(const Malaria_state &y, const double t);
	void update_mosq(const Malaria_state &y, Malaria_state &dydt, const double t);
	void update_human(const Malaria_state &y, Malaria_state &dydt, const double t);
	void update_human(const Human_state &y, Human_state &dydt, const int j, const int k, const double t, const double eird_max);

	void ageing(const Malaria_state &y, Malaria_state &dydt);
	int find_age_group(const double age1) const;
	double mean_in_age_range(const vector<double> &mu, const double age1, const double age2) const;

	void find_equilibrium(Malaria_state &y);	
	void find_equilibrium_known_eir(Malaria_state &y, const double eiry0);

	//******************************************************************************
	//******************************************************************************

	double par(const Human_state &y, const int i, const bool pcr) const;
	void par_prev(const Malaria_state &y, vector<double> &pred, const bool pcr=false) const;
	double par_prev(const Malaria_state &y, const int i, const bool pcr=false) const;
	double par_prev(const Malaria_state &y, const bool pcr=false) const;
	double par_prev(const Malaria_state &y, const double age1, const double age2, const bool pcr=false) const;
	void update_sum_H(const Malaria_state &y);
protected:
	public:
	double increase_IB(const double eir1);
	double increase_IC(const double FOIh1);
	double increase_IV(const double FOIh1);
	double increase_ID(const double FOIh1);
	inline double find_phi(const int i, const int j, const int k, const double ICA);
	inline double find_theta(const int i, const int j, const int k, const double IV);
	inline double find_foi(const double eird_max, const int i, const double IB);
	inline double find_b(const int i, const double IB);
	inline double prob_det(const int i, const double ID) const;

	void find_mean_c(const Malaria_state &y);

	void ageing_and_death(const double *y, double *dydt) const;
	void ageing_mean(const double *y, double *dydt, const double y_at_birth=0.0) const;
	void init_demography(const vector<vector<double> > &v, const double birth_rate0);
	void init_mosq();
	void init_infection();

	//******************************************************************************
	//******************************************************************************
	void human_equilibrium(Human_state &y, const double eird_max, const int j);
	void mosq_equilibrium(Malaria_state &y);
	void update_from_equilibrium(Malaria_state &y);
	void equilibrium_foi(Malaria_state &y_mean, const double eird_max);

	//******************************************************************************
	//******************************************************************************
};

//******************************************************************************
//******************************************************************************

