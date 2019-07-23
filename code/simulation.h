
//******************************************************************************
//******************************************************************************

void change_site_file(const string &file){
	vector<string> s;
	vector<double> x;
	import_file(file, s, x);
	for(int i=0; i<s.size(); i++)
		x[i]=from_map(s[i], -1E20);
	export_file(file, s, x);
}

struct Simulation{
	ofstream* save_out;
	ifstream* save_in;
	void save_state();
	void use_state();

	double runin;
	double final_run;
	string drug_root;
	string demog_file;
	string pop_file;

	Event_manager event_manager;
	Village village;
	vector<Event*> events;
	vector<Repeated_event*> repeated_events;
	int output_run;
	int num_runs;
	int itn_freq;
	int mosq_update_dy;

	Ran sim_ran;

	Village_parms village_parms;
	Human_parms human_parms;
	
	vector<double> epi_ages;
	vector<vector<double> >* to_output;
	vector<Output_var> output_vars;
		
	Simulation() : to_output(0) { }
	~Simulation();
	void setup(const unsigned long long seed0, const int run, vector<vector<double> >* to_output0, vector<Output_var> &output_vars0, const string &root, const string &demog_file0, const string &pop_file0);
	void clear();
	void init(const bool baseline);
	void init_baseline();
	void init_det();
	void init_indiv();

	void add_main_interventions();

	void run();
	double run_det(const double r, const int prev_yrs=0);
	double find_prev();
	double find_time_from_offset(const double offset, const bool offset_absolute);

	bool itn_now();

};

//******************************************************************************
//******************************************************************************
struct Prev_func{
	vector<Simulation>& simulations;
	double prev;
	int prev_yrs;
	Prev_func(vector<Simulation> &s, const double p) : simulations(s), prev(p)  {
		prev_yrs=from_map_int("prev_years", 0, 10, 0);
	}
	double find_prev(const double r){
		double prev1=0;
		const int n=simulations.size();
		const int num_people=from_map_int("num_people");
		const bool feedback=from_map_bool("larval_feedback");
		for(int i=0; i<n; i++){
			Simulation& simulation=simulations[i];
			simulation.village.setup(&simulation, num_people, simulation.demog_file, simulation.pop_file);
		}

		#pragma omp parallel for schedule(dynamic) reduction(+ : prev1)
		for(int i=0; i<n; i++){
			Simulation& simulation = simulations[i];
			if(feedback){
				simulation.village.parms.run_type=2;
				simulation.run_det(r, prev_yrs);
				simulation.village.parms.run_type=1;
			}
			else
				simulation.village.parms.run_type=0;
			prev1 += simulation.run_det(r, prev_yrs);
			for(int m=0; m<num_mosq; m++){
				simulation.village.mosq_pop[m].b_rate_t.clear();
				simulation.village.mosq_pop[m].b_rate.clear();
				simulation.village.mosq_pop[m].birth_index=0;
			}
		}
		prev1 /= n;
		return prev1;
	}
	double operator()(const double r){
		const double prev1=find_prev(exp(r));
		return (sqrt(prev1)-sqrt(prev))/sqrt(max(prev, 0.0001));
	}
};

//******************************************************************************
//******************************************************************************

struct Itn_pulse_func{
	void operator()(Human* h){
		h->give_itn();
	}
};

struct Irs_pulse_func{
	void operator()(Human* h){
		h->give_irs();
	}
};
struct Mda_pulse_func{
	int trt;
	bool screen;
	Mda_pulse_func(const int tr, const bool scr) : trt(tr), screen(scr) {} 
	void operator()(Human* h){
		Village* village=h->village_het->village;
		if(screen)
			village->num_mda_screen++;
		if(!screen || h->uniform()<h->prob_positive(h->find_age(), false)){
			h->treat(trt);
			village->num_mda_trt++;
		}
	}
};
struct Smc_pulse_func{
	int trt;
	Smc_pulse_func(const int tr) : trt(tr) {} 
	void operator()(Human* h){
		h->treat(trt);
		h->village_het->village->num_smc++;
	}
};
struct Pev_pulse_func{
	int vaccine;
	bool boost;
	double boost_cov;
	Pev_pulse_func(const int vacc=2, const bool b=false, const double bc=1.0) : vaccine(vacc), boost(b), boost_cov(bc) {} 
	void operator()(Human* h){
		if(!boost || h->uniform()<boost_cov)
			h->give_pev(vaccine, boost);
	}
};
struct Tbv_pulse_func{
	void operator()(Human* h){
		h->give_tbv();
	}
};

//******************************************************************************
//******************************************************************************

struct Larvicide{
	int num;
	Larvicide(const int n) : num(n) {}
	void operator()(Village* village, const double p) {
		for (int m = 0; m<village->mosq_pop.size(); m++)
			village->mosq_pop[m].change_larvicide(p, num);
	}
};


struct Mass_itn{
	bool check_itn_now;
	int target;
	double age0;
	double age1;
	Mass_itn(const bool check=true, const int tg=0,	const double a0=0.0*dy, const double a1=199.0*dy) : 
		check_itn_now(check), target(tg), age0(a0), age1(a1) {}
	void operator()(Village* village, const double prob){
		if(prob>0 && (!check_itn_now || village->simulation->itn_now()))
			village->det ? village->parms.add_itn(prob) 
				: village->pulse(prob, target, 0, Itn_pulse_func(), age0, age1, true);
	}
};

struct Mass_irs{
	void operator()(Village* village, const double prob){
		if(village->det && prob>0)
			village->parms.add_irs(prob);
		else
			village->pulse(prob, 0, 1, Irs_pulse_func(), 0.0*dy, 200.0*dy, true);
	}
};

struct Mass_smc{
	double age0;
	double age1;
	int trt;
	Mass_smc(const double a0, const double a1, const int trt0) : age0(a0), age1(a1), trt(trt0) {}
	void operator()(Village* village, const double prob){
		village->pulse(prob, 0, 2, Smc_pulse_func(trt), age0, age1, false);
	}
};
//******************************************************************************

struct Change_act_cov{
	void operator()(Village* village, const double prob){
		set_to_zero(village->drug_cov);
		village->drug_cov[1]=prob;
	}
};

struct Change_drug{
	int n;
	Change_drug(const int n0) : n(n0) { }
	void operator()(Village* village, const double prob){
		village->change_drug(n);
	}
};
struct Change_itn{
	int n;
	Change_itn(const int n0) : n(n0) { }
	void operator()(Village* village, const double prob){
		for(int m=0; m<village->mosq_pop.size(); m++)
			village->mosq_pop[m].change_itn(n);
	}
};
struct Change_irs{
	int n;
	Change_irs(const int n0) : n(n0) { }
	void operator()(Village* village, const double prob){
		for(int m=0; m<village->mosq_pop.size(); m++)
			village->mosq_pop[m].change_irs(n);
	}
};
struct Introduce_pev_epi{
	void operator()(Village* village, const double prob){
		village->prob_pev_epi=prob;
	}
};
struct Introduce_ipt_epi{
	int drug;
	Introduce_ipt_epi(const int dr) : drug(dr) { }
	void operator()(Village* village, const double prob){
		village->prob_ipt_epi=prob;
		village->ipt_drug=drug;		
	}
};
struct Introduce_itn_birth{
	double prob_adult;
	Introduce_itn_birth(const double pa) : prob_adult(pa) { } 
	void operator()(Village* village, const double prob){
		if(prob>0){
			village->prob_itn_birth=prob;
			village->prob_itn_birth_adult=prob_adult;
			village->any_itn_irs=true;
		}
	}
};

//******************************************************************************
//******************************************************************************
