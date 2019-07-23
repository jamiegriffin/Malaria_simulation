 //**************************************************************************
//***************************************************************************

struct Village_parms{
	double sigma2;
	vector<double> het_wt;
	vector<double> het_exp;
	double update_interval;
	void setup(Parms &parms);
};

struct Human{
	
	void save_state(ofstream &out);
	void use_state(ifstream &in);
	void save_event(ofstream &out, Event_manager* em, Event &ev);
	void use_event(ifstream &in, Event_manager* em, Event &ev, const double t0);

	Move_state_event move_event;
	Death_event death;
	Leave_itn_event leave_itn_event;
	Epi_event epi_event;
	
	enum Infection_state {S, T, D, A1, U, P};
	
	double date_of_birth;

	double IC_M;
	double IC;
	double IC_time;
	double IC_boost_time;

	double IB;
	double IB_time;
	double IB_boost_time;
	
	double ID;
	double ID_time;
	double ID_boost_time;

	double IV_M;
	double IV;
	double IV_time;
	double IV_boost_time;

	Village_het* village_het;
	Human_parms* hp;
	vector<double> FOIM_contribution;
	
	//***********************************************************************
	// interventions
	ITN_IRS itn_irs;

	//***********************************************************************
	
	bool pev1;
	bool pev2;

	double pev_ab0;
	double pev_d1;
	double pev_d2;
	double pev_rho;

	double pev_time;
	double pev_efficacy;
	
	int num_vacc_doses;
	
	//***********************************************************************
	double tbv_time;
	bool tbv;

	bool ipt;

	vector<double> muj_repeated;
	
	//***********************************************************************
	Infection_state infection_state;

	int trt;
	double rel_c;

	double trt_time;

	//***********************************************************************

	double uniform() const;
	double time_now() const;

	Human(Village_het* vh);
	void set_state(Parms &parms, Human_state &state, const bool itn);

	double find_age() const;
	double find_age(const double t) const;
	int find_age_group(const double age) const;
	
	double prob_infection(const double age);
	double prob_clinical_disease(const double age);
	double prob_severe_disease(const double age);

	double prob_detection(const double age);
	//***********************************************************************
	
	void move_state();
	void handle_new_state(const double t);
	void possible_inf_bite(const double prob_bite, const double t);
	void die();

	void update_FOIM(const bool new_state, const double t);

	double infectivity(const double age);
	double prob_positive(const double age, const bool pcr);

	//***********************************************************************
	
	void give_itn();
	void give_irs();
	void leave_itn();

	void give_epi(const int num_epi);
	void give_pev(const int vaccine, const bool boost);
	void give_tbv();

	void treat(const int new_trt);

};

//***************************************************************************
//***************************************************************************

struct Pev{
	int v;
	void setup(const int v0);
	double find_parameter(const string &p, const double low=0, const double high=1E20, const double x=-99999);

	double find_efficacy(Human& h, const double t);

	bool ab_model;

	// vaccine efficacy parameters as function of anti-body level
	double alpha;
	double beta;
	double Vmax;

	// anti-body decay parameters
	double ab0_mu;
	double ab0_sigma;
	double ab0_boost_mu;
	double ab0_boost_sigma;
	
	double d1_mu;
	double d1_sigma;
	double d2_mu;
	double d2_sigma;
	double rho_mu;
	double rho_sigma;

	double rho_boost_mu;
	double rho_boost_sigma;
	
	double efficacy;
	double decay_scale;
	double decay_scale2;
	double rho;
	double decay_shape;
	double leakiness;

	double boost_efficacy;

	vector<double> efficacy_table;
};

struct Human_parms{

	// model structure choices
	int inf_immunity;
	int det_immunity;
	int clin_immunity;
	int sev_immunity;

	bool separate_sev;
	
	//********************************************************************
	//********************************************************************

	vector<Pev> pevs;
	double pev_epi_boost_coverage;
	
	//********************************************************************
	//********************************************************************

	double tbv_eff;
	double tbv_decay;

	bool leave_itn_on;
	double leave_itn_time;
	bool retain_itn_irs;

	vector<vector<double> > chol_repeated;
	vector<double> sigma_repeated;
	vector<double> temp_prev;

	double eta; // birth and death rate

	double a0;	// relative biting rate at age a is 1-rho*exp(a/a0)
	double rho;
	double br0;

	double meanE; //mean duration in each infection state
	double meanD;
	double meanA;
	double meanU;

	vector<Human::Infection_state> next_state;

	// clinical immunity parameters
	double phi0;
	double phi1;
	double dc;		// duration parameter for acquired clinical immunity
	double kc;
	double IC0;
	double uc;
	double gammac;
	double fc0;
	double ac0;
	double IC00;
	double rc;

	double P_IC_M;	// maternal clinical immunity as proportion of mother's
	double dm;		// duration parameter for maternal clinical immunity

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
	double IB00;
	double rb;

	// detection immunity
	double dmin;
	double dd;
	double kd;
	double ID0;
	double ud;
	double gammad;
	double fd0;
	double ad0;
	double ID00;
	double rd;

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
	double IV00;
	double rv;

	// Infection parameters
	// human to mosquito
	// probability that mosquito becomes infected
	double cD;	
	double cU;	
	double gamma_inf;
	double alpha_pcr;
	double beta_pcr;

	vector<double> fD;
	vector<double> fV;

	void setup(Parms& parms);

	bool use_table;
	static const int Tn=10000;
	vector<double> rel_biting_rate_table;
	vector<double> prob_inf_table;
	vector<double> prob_clin_table;
	vector<double> prob_sev_table;
	vector<double> prob_det_table;

	vector<double> cA_table;
	vector<double> pcrU_table;
	vector<double> pcrA_table;

	double rel_biting_rate(const double age);

	double prob_infection(const double I);
	double find_pcr_par(const double q);
	double prob_clinical_disease(const double I);
	double prob_detection(const double I, const int ai);
	double prob_severe_disease(const double I, const int ai);

	double infectivity(const Human::Infection_state infection_state, const double prob_det, const double rel_c);
	double prob_positive(const Human::Infection_state infection_state, const double prob_det, const bool pcr);
};

//***************************************************************************
//***************************************************************************
struct Village{

	void save_state(ofstream &out);
	void use_state(ifstream &in);

	Human* random_human();
	Guide_table human_guide;

	double num_itn;
	double num_irs;
	double num_trt;
	double num_smc;
	double num_mda_trt;
	double num_mda_screen;
	double num_vaccinees;
	double num_vacc_doses;
	double num_vaccinees_boost;

	bool demog;
	double num_births;
	int year0;

	bool det;
	double peak_time;

	Village_parms* village_parms();
	Parms parms;

	vector<Mosq_pop> mosq_pop;
	vector<Village_het> village_het;
	vector<Drug> drugs;
	vector<double> drug_cov;
	void change_drug(const int n);

	vector<double> sum_repeat;
	vector<double> sum_fail;
	vector<double> sum_repeat_j;
	vector<double> sum_fail_j;

	//***********************************************************************
	Village();
	~Village();
	void clear();

	int N;
	unsigned int num_updates;
	Simulation* simulation;
	
	double EIR_count;
	vector<double> clin_inc;
	vector<double> clin_inc_age0;
	vector<double> clin_inc_age1;

	vector<double> inf_inc;
	vector<double> inf_inc_age0;
	vector<double> inf_inc_age1;

	vector<double> sev_inc;
	vector<double> sev_inc_age0;
	vector<double> sev_inc_age1;

	void update_inc(const double age, const double prob_clin, const double prob_sev);

	void setup(Simulation* sim, const unsigned int N, const string &demog_file, const string &pop_file);
	void init_indiv();

	void update_mosq(const double dt);
	void update_foim();
	void output(const int index, const double time, const double output_interval);
	void mass_birth();

	int total_pop() const;

	//***********************************************************************

	template<class F>
	void pulse(const double prob, const int target, const int j_repeated, F f, const double age0, const double age1, const bool to_update_foim);

	bool any_itn_irs;
	bool any_tbv;

	double prob_pev_epi;
	double prob_ipt_epi;
	double prob_itn_birth;
	double prob_itn_birth_adult;

	int ipt_drug;
	vector<double> old_FOIM_contribution;
	
	//***********************************************************************

	int init_demog_year;
	Guide_table demog_age_guide;
	vector<double> beta;
	vector<vector<double> > alpha;
	vector<double> exit_prob;

};

//***************************************************************************
//***************************************************************************

struct Village_het : public vector<Human*>{
	int het_j;
	Village* village;
	Village_het();
	void setup(Village* v, const int j);
};

//***************************************************************************
//***************************************************************************
