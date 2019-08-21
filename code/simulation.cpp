//***************************************************************
//***************************************************************

int find_map_years(const string &prefix){
	int y=0;
	while(in_map(prefix + to_string(y+1)))
		y++;
	return y;
}
void modify_itn_cov(){

	const int n=find_map_years("itn_cov");
	const double last_itn_cov=n>0 ? from_map("itn_cov"+to_string(n), 0.0, 1.0) : 0.0;
	if(from_map_bool("continue_itn") && !from_map_bool("itn")){
		modify_map("itn_coverage", last_itn_cov);
		modify_map("itn", 1.0);
		modify_map("itn_start", 0.0);
	}		
	
	const int irs_years=find_map_years("irs_cov");
	const double last_irs_cov=irs_years>0 ? from_map("irs_cov"+to_string(irs_years), 0.0, 1.0) : 0.0;
	if(from_map_bool("continue_irs") && !from_map_bool("irs")){
		modify_map("irs_coverage", last_irs_cov);
		modify_map("irs", 1.0);
		modify_map("irs_start", 0.0);
	}

	const bool itn_usage=from_map_bool("itn_usage", false);
	if(itn_usage){
		const int m=from_map_int("num_runs");
		const double D=from_map("itn_leave_dur");
		const double T=from_map("itn_frequency");
		const bool leave_itn=from_map_bool("itn_leave");
		const double r=leave_itn ? exp(-1.0/D) : 1.0;

		vector<double> x(n);
		for(int i=0; i<n; i++){
			const double p0=from_map("itn_cov"+to_string(i+1), 0.0, 1.0);
			double y=0;
			double s=1.0;
			for(int j=i-1; j>=max(i-m+1, 0); j--)
				y += x[j]*(s *= r);
			
			x[i]=m*p0/sqrt(r) - y;
			x[i]=max(0.0, min(1.0, x[i]));
		}
		for(int i=0; i<n; i++)
			modify_map("itn_cov"+to_string(i+1), x[i]);

		if(from_map_bool("itn") && leave_itn){
			const double p=from_map("itn_coverage", 0.0, 1.0);
			const double a=D/T*(1.0-exp(-T/D));
			if(p>a)
				error_crit("Can't have itn_coverage greater than " + to_string(a) + " when it represents usage");
			modify_map("itn_coverage", p/a);
		}
	}
}

void Simulation::setup(const unsigned long long seed0, const int run, vector<vector<double> >* to_output0, 
	vector<Output_var> &output_vars0, const string &root, const string &demog_file0, const string &pop_file0){

	output_vars = output_vars0;
	demog_file=demog_file0;
	pop_file=pop_file0;
	drug_root=root;
	sim_ran.setup(seed0);

	output_run=run;
	num_runs=from_map_int("num_runs");
	itn_freq=from_map_int("itn_frequency", 1);

	to_output=to_output0;
	mosq_update_dy=from_map_int("mosq_update");

	int e=0;
	while(in_map("epi_age_" + to_string(e+1))){
		epi_ages.push_back(from_map("epi_age_" + to_string(e+1), (e==0 ? 0 : epi_ages.back())));
		e++;
	}
	for(int i=0; i<epi_ages.size(); i++)
		epi_ages[i] *= 365.0/12.0;

	const int num_people=from_map_int("num_people");
	village.setup(this, num_people, demog_file, pop_file);
}
double Simulation::find_prev(){
	return village.parms.lag_prev.mean(village.parms.ode.time()-dy);
}

bool Simulation::itn_now(){
	const int y=floor(event_manager.time_now()/dy);
	return modulo(y, itn_freq) == modulo(output_run, itn_freq);
}

void Simulation::init_baseline(){
	const double eps=1E-6;
	const int itn_years=find_map_years("itn_cov");
	const double irs_u=find_time_from_offset(from_map("irs_offset", -1, 1), from_map_bool("irs_offset_absolute"));
	if(itn_years>0){
		if(itn_years*dy > runin)
			error_crit("Number of baseline ITN years must be less than run-in");
		for(int i=0; i<itn_years; i++)
			single_intervention_event(&village, Mass_itn(), from_map("itn_cov"+to_string(i+1), 0, 1), 3*eps + runin -(itn_years-i)*dy);
	}

	const int irs_years=find_map_years("irs_cov");
	if(irs_years>0){
		if(irs_years*dy > runin)
			error_crit("Number of baseline IRS years must be less than run-in");

		for(int i=0; i<irs_years; i++)
			single_intervention_event(&village, Mass_irs(), from_map("irs_cov"+to_string(i+1), 0, 1), 3*eps + runin - (irs_years-i-irs_u)*dy);
	}
}
void Simulation::init_indiv(){
	village.init_indiv();
	repeated_events.push_back(new Mosq_update_event(&village, 1.0/mosq_update_dy));
	repeated_events.push_back(new FOIM_update_event(&village));
	if(village.demog)
		repeated_events.push_back(new Mass_birth_event(&village));
}
void Simulation::init_det(){
	const double h0=from_map("det_step", 0.01, 5.0, 1.0);
	village.parms.ode.setup(village.parms.state, village.parms, -runin, (village.parms.run_type==1 ? h0 : min(h0*0.1, 0.1)));
	village.parms.find_equilibrium(village.parms.state);
}
void Simulation::init(const bool baseline){
	const int num_people=village.N;
	const double span=1*dy;									
	const double width=span/num_people;	
	event_manager.setup(-runin, span, width);
	init_det();
	if(baseline)
		init_baseline();
}
double Simulation::run_det(const double r, const int prev_yrs){
	clear();
	for(int m=0; m<village.mosq_pop.size(); m++)
		village.mosq_pop[m].rel_M=r;

	const int itn_years=find_map_years("itn_cov");
	const int irs_years=find_map_years("irs_cov");
	runin=(village.parms.run_type==2 ? dy*max(15, itn_years+1, irs_years+1) : dy*from_map_int("runin", 10));
	const bool baseline_int=from_map_bool("baseline_int");
	init(baseline_int);

	const double prev_time_dy=dy*prev_yrs;
	event_manager.run_det(runin-prev_time_dy, village.parms.ode);
	return village.parms.run_type<2 ? find_prev() : 0;
}

void Simulation::save_state(){
	sim_ran.save_state(*save_out);
	village.save_state(*save_out);
}
void Simulation::use_state(){
	sim_ran.use_state(*save_in);
	village.use_state(*save_in);
}
void Simulation::run(){

	const int to_save_state=from_map_int("save_state", 0, 2, 0);

	//*******************************************************************************
	const int final_run_yrs=from_map_int("final_run");
	final_run=final_run_yrs*dy;
	const int is=from_map_int("init_output", 0, 20) + 1;
	const double init_store=dy*is;
	const int ri=from_map_int("runin", is + 1);
	runin=dy*ri;
	const double init_ind=max(from_map_int("init_indiv", 0, ri)*dy, init_store + dy);
	const int num_per_yr=from_map_int("output_per_yr", 2);
	const double output_interval=dy/num_per_yr;
	const double eps=1E-6;
	clear();
	
	init(from_map_bool("baseline_int", true));
	repeated_events.push_back(new Output_event(&village, output_interval, runin - init_store - output_interval+2.0*eps));

	if(to_save_state<2){
		if(from_map_bool("det_baseline", false)){
			event_manager.run_det(runin + eps, village.parms.ode);
			return;
		}
		event_manager.run_det(runin - init_ind, village.parms.ode);
		init_indiv();
		event_manager.run(init_ind + eps);
		if(to_save_state==1){
			#pragma omp ordered
			{
				save_state();
			}
		}
	}
	else if(to_save_state==2){
		event_manager.jump_by(runin - init_ind);
		init_indiv();
		event_manager.jump_by(init_ind + eps);
		const double t0=event_manager.time_now();
		for(int i=0; i<events.size(); i++)
			delete events[i];
		events.clear();
		event_manager.clear();	
		event_manager.setup(t0);
		for(int i=0; i<repeated_events.size(); i++)
			repeated_events[i]->skip_to(t0);
		#pragma omp ordered
		{
			use_state();
		}
	}
	event_manager.run(3.0*eps);

	// it's now time 0
	add_main_interventions();

	event_manager.run(final_run-eps);
}

double Simulation::find_time_from_offset(const double offset, const bool offset_absolute){
	double u=offset;
	if(!offset_absolute)
		u += village.peak_time;
	while(u<0)
		u += 1;
	while(u>1)
		u -= 1;
	return u;
}
void Simulation::add_main_interventions(){
	village.num_itn=0;
	village.num_irs=0;
	village.num_trt=0;
	village.num_smc=0;
	village.num_mda_trt=0;
	village.num_mda_screen=0;
	village.num_vacc_doses=0;
					 
	const bool itn=from_map_bool("itn", false);
	const bool irs=from_map_bool("irs", false);
	const bool mda=from_map_bool("mda", false);
	const bool mass_pev=from_map_bool("mass_pev", false);
	const bool epi_pev=from_map_bool("epi_pev", false);
	const bool epi_ipt=from_map_bool("epi_ipt", false);
	const bool smc=from_map_bool("smc", false);
	const bool change_drug=from_map_bool("change_drug", false);
	const bool tbv=from_map_bool("tbv", false);
	const bool itn_birth=from_map_bool("itn_birth", false);

	const bool larvicide = from_map_bool("larvicide", false);

	const bool change_itn=from_map_bool("change_itn", false);
	const bool change_irs=from_map_bool("change_irs", false);

	if(change_drug){
		const string prefix="change_drug";
		string suffix="";
		int num=1;
		do{
			single_intervention_event(&village, Change_drug(num), 0, from_map(prefix + "_time" + suffix)*dy); 
			suffix="_" + to_string(++num);
		}while(in_map(prefix + suffix) && from_map_bool(prefix + suffix));
	}

	if(change_itn){
		const string prefix="change_itn";
		string suffix="";
		int num=1;
		do{
			single_intervention_event(&village, Change_itn(num), 0, from_map(prefix + "_time" + suffix)*dy); 
			suffix="_" + to_string(++num);
		}while(in_map(prefix + suffix) && from_map_bool(prefix + suffix));
	}
	if(change_irs){
		const string prefix="change_irs";
		string suffix="";
		int num=1;
		do{
			single_intervention_event(&village, Change_irs(num), 0, from_map(prefix + "_time" + suffix)*dy); 
			suffix="_" + to_string(++num);
		}while(in_map(prefix + suffix) && from_map_bool(prefix + suffix));
	}	

	if(itn){
		const string prefix="itn";
		repeated_intervention_event(&village, Itn_pulse_func(), prefix, 0, 0, true); 
	}

	if(!itn && from_map_bool("itn_flexible", false)){
		const int ny=from_map_int("final_run");
		for(int y=0; y<ny; y++){
			const string suffix="_" + to_string(y) + "_" + to_string(output_run);
			const double cov=from_map("itn_cov" + suffix, 0, 1);
			const int target=from_map_int_suffix("itn_target", suffix, 0, 3, 0);
			const double age0=from_map_suffix("itn_age0", suffix, 0.0, 200.0, 0.0)*dy;
			const double age1=from_map_suffix("itn_age1", suffix, 0.0, 200.0, 199.0)*dy;
			if(cov>0)
				single_intervention_event(&village, Mass_itn(false, target, age0, age1), cov, y*dy);
		}
	}

	if(itn_birth){
		single_intervention_event(&village, Introduce_itn_birth(from_map("itn_birth_cov_adult", 0, 1)), from_map("itn_birth_coverage", 0, 1), 0.0);
	}
	if(irs){
		const string prefix="irs";
		const bool offset_absolute=from_map_bool(prefix + "_offset_absolute");
		string suffix="";
		int num=1;
		do{
			const double u=find_time_from_offset(from_map_suffix(prefix + "_offset", suffix, -1, 1), offset_absolute);
			repeated_intervention_event(&village, Irs_pulse_func(), prefix, 1, u*dy, true, suffix); 
			suffix="_" + to_string(++num);
		}while(in_map(prefix + suffix) && from_map_bool(prefix + suffix));
	}

	if(mda){
		const string prefix="mda";
		const bool offset_absolute=from_map_bool(prefix + "_offset_absolute");
		string suffix="";
		int num=1;
		do{
			const bool screen=from_map_bool_suffix(prefix + "_screen", suffix);
			const int drug=from_map_int_suffix(prefix + "_drug", suffix);
			const double u=find_time_from_offset(from_map_suffix(prefix + "_offset", suffix, -1, 1), offset_absolute);
			repeated_intervention_event(&village, Mda_pulse_func(drug, screen), prefix, 2, u*dy, false, suffix); 
			suffix="_" + to_string(++num);
		}while(in_map(prefix + suffix) && from_map_bool(prefix + suffix));
	}

	if(smc){
		const string prefix="smc";
		string suffix="";
		int num=1;
		const int num_smc=from_map_int("num_smc", 1, 12, 3); 
		const double smc_offset=from_map("smc_offset", -1.0, 1.0, 0.0);
		do{
			for(int i=0; i<num_smc; i++){
				double u=village.peak_time + (i-(num_smc-1.0)/2.0)/12.0 + smc_offset;
				while(u<=0)
					u += 1;
				const int drug=from_map_int_suffix(prefix + "_drug", suffix);
				repeated_intervention_event(&village, Smc_pulse_func(drug), prefix, 2, u*dy, false, suffix); 
			}
			suffix="_" + to_string(++num);
		}while(in_map(prefix + suffix) && from_map_bool(prefix + suffix));
	}
	
	
	if(mass_pev){
		const string prefix="mass_pev";
		string suffix="";
		int num=1;
		
		const double dose_interval=from_map("mass_pev_dose_interval", 0, 1E20, 1)*dy/12.0;
		vector<double> dose_intervals = { 0.0, dose_interval, 2.0*dose_interval };
		if (from_map_bool("mass_pev_boost1", false))
			dose_intervals.push_back(2.0*dose_interval + from_map("mass_pev_boost1_interval", 0, 1E20, 12)*dy/12.0);
		if (from_map_bool("mass_pev_boost2", false))
			dose_intervals.push_back(2.0*dose_interval + from_map("mass_pev_boost2_interval", 0, 1E20, 24)*dy/12.0);
		sort(dose_intervals.begin(), dose_intervals.end());
		const double boost_coverage = from_map("mass_pev_boost_coverage", 0, 1, 0.0);
		do{
			repeated_intervention_event(&village, Pev_pulse_func(2, boost_coverage, dose_intervals), prefix, 3, 0.0, false, suffix);
			suffix="_" + to_string(++num);
		}while(in_map(prefix + suffix) && from_map_bool(prefix + suffix));
	}
	if(epi_pev){
		const string prefix="pev_epi";
		const int pev_epi_scale_up=from_map_int("pev_epi_scale_up");
		if(pev_epi_scale_up==0)
			single_intervention_event(&village, Introduce_pev_epi(), from_map("pev_epi_coverage", 0, 1), from_map("pev_epi_start")*dy);
		else
			for(int y=0; y<pev_epi_scale_up; y++)
				single_intervention_event(&village, Introduce_pev_epi(), from_map("pev_epi_coverage" + to_string(y+1), 0, 1), (from_map("pev_epi_start")+y)*dy);
	}
	if(epi_ipt){
		single_intervention_event(&village, Introduce_ipt_epi(from_map_int("ipt_drug")), from_map("ipt_epi_coverage", 0, 1), from_map("ipt_epi_start")*dy);
	}
	if(tbv){
		const string prefix="tbv";
		string suffix="";
		int num=1;
		do{
			repeated_intervention_event(&village, Tbv_pulse_func(), prefix, 3, 0, true, suffix); 
			suffix="_" + to_string(++num);
		}while(in_map(prefix + suffix) && from_map_bool(prefix + suffix));
	}

	if(larvicide){
		const string prefix = "larvicide";
		string suffix = "";
		int num = 1;
		do {
			const double larvicide_on = from_map_suffix("larvicide_on", suffix)*dy;
			const double larval_interval = from_map_suffix("larval_interval", suffix)*dy;
			single_intervention_event(&village, Larvicide(num), 1.0, larvicide_on);
			single_intervention_event(&village, Larvicide(num), 0.0, larvicide_on + larval_interval);
			suffix = "_" + to_string(++num);
		} while (in_map(prefix + suffix) && from_map_bool(prefix + suffix));
	}

}

void Simulation::clear(){
	event_manager.clear();
	village.clear();
	for(int i=0; i<events.size(); i++)
		delete events[i];
	for(int i=0; i<repeated_events.size(); i++)
		delete repeated_events[i];
	events.clear();
	repeated_events.clear();
}

Simulation::~Simulation(){
	for(int i=0; i<events.size(); i++)
		delete events[i];
	for(int i=0; i<repeated_events.size(); i++)
		delete repeated_events[i];
	events.clear();
	repeated_events.clear();
}

//***************************************************************
//***************************************************************

double recalculate_M(const bool mult_runs, 
				const string &root, const string &demog_file, const string &pop_file){
	// calculate M from baseline prevalence, including possible ITNs

	const double prev0=from_map("prev", 0, 1);

	const int num_runs=mult_runs ? from_map_int("num_runs") : 1;
	vector<Simulation> simulations(num_runs); 
	vector<Output_var> v(0);
	for(int run=0; run<num_runs; run++)
		simulations[run].setup(1, run, 0, v, root, demog_file, pop_file);
	
	const double total_M_input = from_map("total_M");
	const double prev_tol=from_map("prev_tol", 1.0E-8, 5.0E-2);
	Prev_func f(simulations, prev0);
	const double r=find_root(f, 0.0, prev_tol);
	return total_M_input*exp(r);
}

void run_final(vector<vector<vector<double> > > &to_output, vector<Output_var> &output_vars, 
			const string &root, const string &demog_file, const string &pop_file, const string &saved_file){
	const bool det_baseline=from_map_bool("det_baseline");
	const int itn_years=find_map_years("itn_cov");
	const int irs_years=find_map_years("irs_cov");
	const int num_runs=itn_years==0 && irs_years==0 && det_baseline ? 1 : to_output.size();
	Ran init_ran(from_map_ll("seed"));

	ofstream save_out;
	ifstream save_in;
	const int to_save_state=from_map_int("save_state", 0, 2, 0);
	if(to_save_state==1){
		save_out.open(saved_file.c_str(), ios::out | ios::binary);
		if(!save_out)
			error_crit("Can't open file " + saved_file + " for saving state");
	}
	else if(to_save_state==2){
		save_in.open(saved_file.c_str(), ios::in | ios::binary);
		if(!save_in)
			error_crit("Can't open file " + saved_file + " for using saved state");
	}
	vector<Simulation> simulations(num_runs); 
	for(int run=0; run<num_runs; run++){
		simulations[run].setup(init_ran.ran64(), run, &to_output[run], output_vars, root, demog_file, pop_file);
		if(to_save_state==1)
			simulations[run].save_out=&save_out;
		if(to_save_state==2)
			simulations[run].save_in=&save_in;
	}
	const bool feedback=from_map_bool("larval_feedback");

	#pragma omp parallel for ordered schedule(static,1)
	for(int run=0; run<num_runs; run++){
		Simulation& simulation=simulations[run];
		Village& v=simulation.village;
		if(feedback){
			v.parms.run_type=2;
			simulation.run_det(1);
			v.parms.run_type=1;
		}
		else
			v.parms.run_type=0;
		simulation.run();	
	}
}

//***************************************************************
//***************************************************************

void output_sim(const vector<vector<vector<double> > > &to_output, const string &output_file, 
				const vector<string> &var_names, const string &name, const int final_output_year){

	const int num_vars=var_names.size();
	const int num_num=to_output[0].size() - num_vars - 1;
	const int num_per_yr=from_map_int("output_per_yr");
	const int to_save_state=from_map_int("save_state", 0, 2, 0);
	vector<vector<double> > to_output_mean(to_output[0].size(), vector<double>(to_output[0][0].size()));
	for(int t=0; t<to_output_mean[0].size(); t++){
		to_output_mean[0][t]=to_output[0][0][t];
		for(int k=1; k<to_output_mean.size(); k++){
			to_output_mean[k][t]=0;
			for(int run=0; run<to_output.size(); run++)
				to_output_mean[k][t] += to_output[run][k][t];
			if(k<to_output_mean.size()-num_num)
				to_output_mean[k][t] /= to_output.size();
		}
	}
	vector<vector<double> > to_output_smooth(to_output_mean.size(), 
			vector<double>(to_output_mean[0].size()-num_per_yr, 0));
	for(int k=1; k<to_output_mean.size(); k++){
		for(int t=0; t<num_per_yr; t++)
			to_output_smooth[k][0] += to_output_mean[k][t]/num_per_yr;
		for(int t1=1; t1<to_output_smooth[k].size(); t1++)
			to_output_smooth[k][t1] = to_output_smooth[k][t1-1] + 
				(to_output_mean[k][t1+num_per_yr-1]-to_output_mean[k][t1-1])/num_per_yr;
	}

	//*******************************************************************************
	const int output_header=from_map_int("output_header", 0, 2);
	const int output_type=from_map_int("output_type", 0, 1);
	const bool output_binary=from_map_bool("output_binary", false);

	ofstream out;
	if(output_binary && output_header<2)
		out.open(output_file.c_str(), ios::out | ios::binary);
	else
		out.open(output_file.c_str());

	Output output(out, output_binary);

	if((output_header==1 && !output_binary) || output_header==2){
		if(output_type>0)
			out << "name" << '\t';
		out << "year";
		for(int i=0; i<var_names.size(); i++){
			out << '\t' << var_names[i];
			if(output_type==0)
				out << '\t' << var_names[i] + "_smooth";
		}
		out << "\tnum_itn\tnum_irs\tnum_trt\tnum_smc\tnum_mda_trt\tnum_mda_screen\tnum_vacc_doses\n";
	}

	if(output_header<2){
		if(output_type==0){
			for(int t=num_per_yr; t<to_output_mean[0].size(); t++){
				const double year=to_output_mean[0][t]/dy;
				if(to_save_state==2 && year<0.0)
					continue;
				output(year);
				for(int k=1; k<to_output_mean.size()-num_num; k++){
					const double x0=to_output_mean[k][t];
					output((x0>1E-10 ? x0 : (x0>=-1E-10 ? 0 : -999)));
					const double x1=(to_save_state==2 && year<1.0) ? -999 : to_output_smooth[k][t-num_per_yr];
					output((x1>1E-10 ? x1 : (x1>=-1E-10 ? 0 : -999)));
				}
				for(int k=to_output_mean.size()-num_num; k<to_output_mean.size(); k++)
					output(to_output_mean[k][t]);
				output.endl();
			}
		}
		else{
			for(int y=0; y<=final_output_year; y++){
				const double year=from_map("output_year" + to_string(y), -20);
				if(to_save_state==2 && year<1.0)
					continue;
				for(int t=num_per_yr; t<to_output_mean[0].size(); t++){
					if(abs(to_output_mean[0][t]/dy - year)<0.0005){	
						output(name);
						output(year);
						for(int k=1; k<to_output_mean.size()-num_num; k++){
							const double x1=to_output_smooth[k][t-num_per_yr];
							output((x1>1E-10 ? x1 : (x1>=-1E-10 ? 0 : -999)));
						}
						for(int k=to_output_mean.size()-num_num; k<to_output_mean.size(); k++)
							output(to_output_mean[k][t]);
						output.endl();
						break;
					}
				}
			}
		}
	}

	out.close();

}
//***************************************************************
//***************************************************************
