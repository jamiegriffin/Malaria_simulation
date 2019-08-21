//***************************************************************************
//***************************************************************************

Village_parms* Village::village_parms(){
	return &simulation->village_parms;
}

//***************************************************************************

void Village_parms::setup(Parms& parms1){
	update_interval=1.0/from_map_int("mosq_update");
	sigma2=parms1.sigma2;
	het_wt=parms1.het_wt;
	het_exp=parms1.rel_foi;
}

double Pev::find_parameter(const string &p, const double low, const double high, const double x){
	const string s=to_string(v);
	return in_map("pev" + s + "_" + p) ? from_map("pev" + s + "_" + p, low, high) : from_map("pev_" + p, low, high, x);
}
void Pev::setup(const int v0){
	v=v0;
	const string s=to_string(v);
	ab_model=in_map("ab_model"+s) ? from_map_bool("ab_model"+s) : from_map_bool("ab_model", false);
	
	const double X=1E20;
	if(ab_model){
		alpha = find_parameter("alpha");
		beta = find_parameter("beta");
		Vmax = find_parameter("Vmax", 0, 1);
 
		ab0_mu = find_parameter("ab_mu", -X);
		d1_mu = find_parameter("d1_mu", -X);
		d2_mu = find_parameter("d2_mu", -X);
		rho_mu = find_parameter("rho_mu", -X);
 		ab0_sigma = find_parameter("ab_sigma");
		d1_sigma = find_parameter("d1_sigma");
		d2_sigma = find_parameter("d2_sigma");
		rho_sigma = find_parameter("rho_sigma");
 
		ab0_boost_mu = find_parameter("ab_boost_mu", -X, X, ab0_mu);
		rho_boost_mu = find_parameter("rho_boost_mu", -X, X, rho_mu);
 
		ab0_boost_sigma = find_parameter("ab_boost_sigma", 0, X, ab0_sigma);
		rho_boost_sigma = find_parameter("rho_boost_sigma", 0, X, rho_sigma);

		const int Tn=10000;
		efficacy_table.resize(Tn);
		for(int i=0; i<Tn; i++)
			efficacy_table[i]=Vmax*(1.0-1.0/(1.0+pow(((i+0.5)/Tn)*100.0, alpha)));
	}
	else{
		efficacy = find_parameter("efficacy", 0, 1);
		decay_shape = find_parameter("decay_shape", 0, 1000, 1.0);
		
		const double half_life=find_parameter("half_life");
		const double half_life2=find_parameter("half_life2", 0, X, 0.0001);
		decay_scale = 1/(dy*half_life)*pow(-log(0.5), 1/decay_shape);
		decay_scale2 = 1/(dy*half_life2)*log(2.0);
		
		rho = find_parameter("rho", 0, 1, 1);
		leakiness = find_parameter("leakiness", 0, 10000, 10000); 
		boost_efficacy = find_parameter("boost_efficacy", 0, 1, efficacy);
		
		efficacy=max(1E-8, min(1.0-1E-8, efficacy));
		boost_efficacy=max(1E-8, min(1.0-1E-8, boost_efficacy));
	}
}
double Pev::find_efficacy(Human& h, const double t_ago){
	if(ab_model){
		const double log2=0.6931471805599453;
		const double ab=h.pev_ab0*(h.pev_rho*exp(-t_ago*log2/h.pev_d1)+(1-h.pev_rho)*exp(-t_ago*log2/h.pev_d2));
		const int Tn=efficacy_table.size();
		return efficacy_table[int(min(ab*Tn/(100.0*beta), Tn-1.0))];
	}
	else{
		const double r1=t_ago*decay_scale;
		const double r2=t_ago*decay_scale2;
		return h.pev_efficacy*(decay_shape==1 ? (rho*exp(-r1) + (1.0-rho)*exp(-r2)) : exp(-pow(r1, decay_shape)));
	}
}

void Human_parms::setup(Parms& parms){
	
	inf_immunity =	parms.inf_immunity;
	clin_immunity = parms.clin_immunity;
	sev_immunity = parms.sev_immunity;
	det_immunity = parms.det_immunity;

	if((inf_immunity!=4 && inf_immunity!=6) || (clin_immunity!=4 && clin_immunity!=6) || 
			(sev_immunity!=4 && sev_immunity!=6) || (det_immunity!=4 && det_immunity!=6)){
		error_crit("inf_immunity, clin_immunity, sev_immunity and det_immunity must be 4 or 6 (other options removed in Jan 2017 to simplify code)");
	}

	separate_sev=parms.separate_sev;

	leave_itn_on=from_map_bool("itn_leave");
	leave_itn_time=dy*from_map("itn_leave_dur");

	retain_itn_irs=true;

	pevs.resize(2);
	pevs[0].setup(1);
	pevs[1].setup(2);

	pev_epi_boost_coverage= from_map("pev_epi_boost_coverage", 0, 1, 0.0);
	
	tbv_eff=from_map("tbv_efficacy", 0, 1);
	tbv_decay=log(2.0)/(dy*from_map("tbv_half_life"));

	eta=parms.eta;
	
	meanE=parms.dur_E;
	a0=parms.a0;
	rho=parms.rho;
	br0 = Tn/(100.0*365.0);

	phi0=parms.phi0;
	phi1=parms.phi1;
	dc=parms.dc;
	kc=parms.kc;
	IC0=parms.IC0;
	uc=parms.uc;
	gammac=parms.gammac;
	fc0=parms.fc0;
	ac0=parms.ac0;
	P_IC_M=parms.P_IC_M;
	dm=parms.dm;
	IC00 = Tn/(200.0*IC0);
	rc = 1.0/dc;

	bh=parms.bh;
	bmin=parms.bmin;
	db=parms.db;
	kb=parms.kb;
	IB0=parms.IB0;
	ub=parms.ub;
	gammab=parms.gammab;
	fb0=parms.fb0;
	ab0=parms.ab0;
	IB00 = Tn/(15.0*IB0);
	rb = 1.0/db;

	dmin=parms.dmin;
	dd=parms.dd;
	kd=parms.kd;
	ID0=parms.ID0;
	ud=parms.ud;
	gammad=parms.gammad;
	fd0=parms.fd0;
	ad0=parms.ad0;
	ID00 = Tn/(200.0*ID0);
	rd = 1.0/dd;

	theta0=parms.theta0;
	theta1=parms.theta1;
	dv=parms.dv;
	kv=parms.kv;
	IV0=parms.IV0;
	uv=parms.uv;
	gammav=parms.gammav;
	fv0=parms.fv0;
	av0=parms.av0;
	P_IV_M=parms.P_IV_M;
	dvm=parms.dvm;
	tau_v=parms.tau_v;
	IV00 = Tn/(50.0*IV0);
	rv = 1.0/dv;

	gamma_inf=parms.gamma_inf;
	alpha_pcr=parms.alpha_pcr;
	beta_pcr=parms.beta_pcr;
	meanA=from_map("dur_A");
	meanD=parms.dur_D;
	meanU=parms.dur_U;
	cU=parms.cU;
	cD=parms.cD;
	fD=parms.fD;
	fV=parms.fV;

	const int num_repeated=5;
	chol_repeated=vector<vector<double> >(num_repeated, vector<double>(num_repeated, 0.0));
	sigma_repeated.resize(num_repeated);
	vector<vector<double> > var_repeated=chol_repeated;

	const string names[num_repeated]={"itn","irs","mda","vacc","epi"};
	for(int i=0; i<num_repeated; i++){
		const string name=names[i];
		double corr=from_map(name+"_corr", 0.0, 1.0, 1.0);
		if(corr==0)
			corr=0.00000001;
		if(corr==1)
			corr=0.9999999;
		sigma_repeated[i]=sqrt(corr/(1-corr));
		var_repeated[i][i]=sigma_repeated[i]*sigma_repeated[i];
	}

	for(int i=0; i<num_repeated; i++){
		for(int j=i+1; j<num_repeated; j++){
			const string si=names[i];
			const string sj=names[j];
			double corr=0.0;
			if(in_map(si +"_" + sj +"_corr"))
				corr=from_map(si +"_" + sj +"_corr", -1.0, 1.0);
			else if(in_map(sj +"_" + si +"_corr"))
				corr=from_map(sj +"_" + si +"_corr", -1.0, 1.0);
			if(corr==1.0 || corr==-1.0)
				corr *= 0.9999999;
			var_repeated[i][j]=sqrt(var_repeated[i][i]*var_repeated[j][j])*corr;
			var_repeated[j][i]=var_repeated[i][j];		
		}
	}

	cholesky(var_repeated, chol_repeated, 1E-12);

	const int num_states=parms.state.human_state[0][0].num_prev_states();
	temp_prev.resize(num_states);

	next_state.resize(num_states);

	next_state[Human::T]=from_map_bool("p_protect") ? Human::S : Human::P;
	next_state[Human::P]=Human::S;
	next_state[Human::D]=Human::A1;
	next_state[Human::A1]=Human::U;
	next_state[Human::U]=Human::S;
	use_table=true;

	prob_inf_table.resize(Tn);
	for(int i=0; i<Tn; i++)
		prob_inf_table[i]=(inf_immunity==0 ? bh : (bh*(1-bmin)/(1+pow((i+0.5)*15.0/Tn, kb)) + bh*bmin));

	prob_clin_table.resize(Tn);
	for(int i=0; i<Tn; i++)
		prob_clin_table[i]=phi0*((1-phi1)/(1+pow((i+0.5)*200.0/Tn, kc)) + phi1);
	
	rel_biting_rate_table.resize(Tn);
	for(int i=0; i<Tn; i++)
		rel_biting_rate_table[i]=1-rho*exp(-(i+0.5)*365.0*100.0/(Tn*a0));
				
	prob_sev_table.resize(Tn);
	for(int i=0; i<Tn; i++)
		prob_sev_table[i]=pow((i+0.5)*50.0/Tn, kv);

	prob_det_table.resize(Tn);
	for(int i=0; i<Tn; i++)
		prob_det_table[i]=pow((i+0.5)*200.0/Tn, kd);
			
	cA_table.resize(Tn);
	for(int i=0; i<Tn; i++)
		cA_table[i]= cU + (cD-cU)*pow((i+0.5)/Tn, gamma_inf);
	
	pcrU_table.resize(Tn);
	for(int i=0; i<Tn; i++)
		pcrU_table[i]=pow((i+0.5)/Tn, beta_pcr);
		
	pcrA_table.resize(Tn);
	for(int i=0; i<Tn; i++)
		pcrA_table[i]=pow((i+0.5)/Tn, alpha_pcr);

}

//***************************************************************************
//***************************************************************************

Village::Village() : EIR_count(0), num_updates(0), any_itn_irs(false), any_tbv(false), prob_pev_epi(0), prob_ipt_epi(0),
		prob_itn_birth(0), prob_itn_birth_adult(0), ipt_drug(-1), det(true), peak_time(0),
		num_itn(0), num_irs(0), num_trt(0), num_smc(0), num_mda_trt(0), num_mda_screen(0), 
		num_vacc_doses(0), num_births(0)  { }

Village::~Village(){
	for(int j=0; j<village_het.size(); j++)
		for(int i=0; i<village_het[j].size(); i++)
			delete village_het[j][i];
	village_het.clear();
}

void Village::clear(){
	for(int j=0; j<village_het.size(); j++)
		for(int i=0; i<village_het[j].size(); i++)
			delete village_het[j][i];
	village_het.clear();
	parms.itn_irs_det.itn=false;
	parms.itn_coverage=0;
	parms.itn_irs_det.irs=false;
	parms.irs_coverage=0;

	for(int m=0; m<mosq_pop.size(); m++)
		mosq_pop[m].reset();
	parms.state.current_size=parms.state.init_size;
	parms.clear_history();
}

void Village::init_indiv(){

	det=false;
	old_FOIM_contribution.resize(mosq_pop.size());
	clin_inc.assign(clin_inc.size(), 0);
	inf_inc.assign(inf_inc.size(), 0);
	sev_inc.assign(sev_inc.size(), 0);

	for(int m=0; m<mosq_pop.size(); m++)
		mosq_pop[m].C=0;

	village_het.resize(village_parms()->het_exp.size());
	for(int j=0; j<village_het.size(); j++)
		village_het[j].setup(this, j);
	const double itn_time=(parms.itn_irs_det.itn ? parms.itn_irs_det.itn_time : 0);
	try{
		vector<double> het_u(N);
		for(int i=0; i<N; i++)
			het_u[i]=(i+uniform(simulation->sim_ran))/N;
		for(int i=1; i<N; i++)
			swap(het_u[i], het_u[ran_int(i, simulation->sim_ran)]);
		for(int i=0; i<N; i++){
			const int j=(i<village_het.size() ? i : discrete_dev(village_parms()->het_wt, het_u[i], true));
			village_het[j].push_back(new Human(&village_het[j]));
		}
	}
	catch(bad_alloc){
		error_crit("Ran out of memory in Village::init_indiv");
	}
	if(parms.itn_coverage>0)
		pulse(parms.itn_coverage, 0, 0, Itn_pulse_func(), 0.0*dy, 200.0*dy, false);

	for(int j=0; j<village_het.size(); j++){
		for(int i=0; i<village_het[j].size(); i++){
			Human* h=village_het[j][i];
			const bool itn=h->itn_irs.itn;
			if(itn)
				h->itn_irs.itn_time=itn_time;
			const int k=(itn ? 1 : 0);
			h->set_state(parms, parms.state.human_state[k][j], itn);
		}
	}
	
	vector<double> g(village_het.size());
	for(int j=0; j<village_het.size(); j++)
		g[j]=village_het[j].size();
	human_guide.setup(g);

	const double t=simulation->event_manager.time_now();
	for(int m=0; m<mosq_pop.size(); m++)
		mosq_pop[m].set_init_indiv(parms.state.mosq_state[m], t);

	update_foim();
}

void Village::save_state(ofstream &out){
	output_binary(num_updates, out);
	output_binary(EIR_count, out);
	output_binary(any_itn_irs, out);
	::save_state(clin_inc, out);
	::save_state(inf_inc, out);
	::save_state(sev_inc, out);
	for(int m=0; m<mosq_pop.size(); m++)
		mosq_pop[m].save_state(out);
	for(int j=0; j<village_het.size(); j++)
		for(int i=0; i<village_het[j].size(); i++)
			village_het[j][i]->save_state(out);
}

void Village::use_state(ifstream &in){
	input_binary(num_updates, in);
	input_binary(EIR_count, in);
	input_binary(any_itn_irs, in);
	::use_state(clin_inc, in);
	::use_state(inf_inc, in);
	::use_state(sev_inc, in);
	for(int m=0; m<mosq_pop.size(); m++)
		mosq_pop[m].use_state(in);
	for(int j=0; j<village_het.size(); j++)
		for(int i=0; i<village_het[j].size(); i++)
			village_het[j][i]->use_state(in);
}
void Village::setup(Simulation* sim, const unsigned int N0, const string &demog_file, const string &pop_file){
	simulation=sim;
	parms.simulation=sim;
	N=N0;

	mosq_pop.assign(num_mosq, this);
	parms.mosq_pop = &mosq_pop;
	const int num_het=from_map_int("num_het");
	const int num_age=from_map_int("num_age");

	sum_repeat.resize(num_mosq);
	sum_fail.resize(num_mosq);
	sum_repeat_j.resize(num_mosq);
	sum_fail_j.resize(num_mosq);

	demog = demog_file.size()>0;
	vector<vector<double> > v(2);
	if(demog){
		alpha.clear();
		beta.clear();
		year0=from_map_int("year0", 1900, 2100, 2014);
		vector<vector<double> > u;
		import_file(demog_file, u);

		init_demog_year=u[0][0];

		for(int i=0; i<u.size(); i++){
			if(u[i][0]==init_demog_year){
				v[0].push_back(u[i][1]);
				v[1].push_back(u[i][4]);
			}
		}

		int previous_year=init_demog_year-1;
		for(int i=0; i<u.size(); i++){
			const int year=u[i][0];
			if(year != previous_year){
				beta.push_back(u[i][3]);
				alpha.push_back(vector<double>());
				previous_year=year;
			}
			alpha.back().push_back(u[i][4]);
		}
		vector<double> w(v[0].size(), 5.0);
		for(int i=0; i<int(w.size())-2; i++)
			w[i]=v[0][i+1]-v[0][i];
		demog_age_guide.setup(w);

		for(int y=0; y<alpha.size(); y++){
			const double ma=max_el(alpha[y]);
			for(int i=0; i<alpha[y].size(); i++)
				alpha[y][i] /= ma;
		}
	}
	else{
		const double eta=from_map("eta")*365.0;
		v[0].assign(1, 0.0);
		v[1].assign(1, eta);
		beta=v[1];
	}

	parms.setup(num_age, num_het, v, beta[0]);

	simulation->village_parms.setup(parms);
	simulation->human_parms.setup(parms);

	const string mosq[num_mosq]={"fun","arab","gamb_ss"};
	for(int m=0; m<mosq_pop.size(); m++)
		mosq_pop[m].setup(&parms, mosq[m], N, pop_file);

	peak_time=0;
	double max_f=-1E20;
	for(int i=0; i<365; i++){
		const double f=mosq_pop[0].seasonal_curve(i+0.5);
		if(f>max_f){
			max_f=f;
			peak_time=(i+0.5)/365.0;
		}
	}
	
	drugs.clear();
	int drug=0;
	string drug_str="drug" + to_string(drug);
	do{
		const string drug_file=simulation->drug_root +"p_protect_"+to_string(drug)+".txt";
		drugs.push_back(Drug());
		drugs.back().setup(from_map("dur_T"), from_map(drug_str +"_rel_c", 0, 1), from_map(drug_str +"_eff", 0, 1), drug_file);
		if(!from_map_bool("p_protect")){
			drugs.back().dur_P=from_map(drug_str +"_dur_P");
			drugs.back().shape_P=from_map(drug_str +"_shape_P", 0.01, 1E20, 1.0);
		}
		drug_str="drug" + to_string(++drug);
	}while(in_map(drug_str +"_eff"));
	drug_cov.assign(drugs.size(), 0);
	change_drug(0);
	parms.rec_T=1.0/from_map("dur_T");
	if(max_el(drug_cov)==0.0){
		parms.ft=0.0;
		parms.rec_P=1.0/max(0.1, drugs[0].dur_P-drugs[0].dur_T);
		parms.cT=drugs[0].rc*parms.cD;
	}
	else{
		parms.rec_P=0;
		parms.cT=0;
		parms.ft=0;
		for(int i=0; i<drugs.size(); i++){
			const double wi=drug_cov[i]*drugs[i].efficacy;
			parms.ft += wi;
			parms.rec_P	+= wi*max(0.1, drugs[i].dur_P-drugs[i].dur_T);
			parms.cT += wi*drugs[i].rc;
		}
		parms.rec_P=parms.ft/parms.rec_P;
		parms.cT *= parms.cD/parms.ft;
	}

	for(int l=0; l<simulation->output_vars.size(); l++){
		Output_var* o=&(simulation->output_vars)[l];
		if(o->quantity==Output_var::clin_inc){
			clin_inc_age0.push_back(o->age0);
			clin_inc_age1.push_back(o->age1);
		}
		if(o->quantity==Output_var::sev_inc){
			sev_inc_age0.push_back(o->age0);
			sev_inc_age1.push_back(o->age1);
		}
		if(o->quantity==Output_var::inf_inc){
			inf_inc_age0.push_back(o->age0);
			inf_inc_age1.push_back(o->age1);
		}
	}
	inf_inc.resize(inf_inc_age0.size());
	clin_inc.resize(clin_inc_age0.size());
	sev_inc.resize(sev_inc_age0.size());
}

void Village::change_drug(const int n){
	const string s="_" + to_string(n);
	if(n==0 && in_map("ft") && in_map("act_baseline_cov") && !in_map("drug_cov_" + to_string(0) + s)){
		set_to_zero(drug_cov);
		const double f=from_map("ft", 0, 1);
		const double q=from_map("act_baseline_cov", 0, 1);
		drug_cov[0]=f*(1-q);
		drug_cov[1]=f*q;
		return;
	}
	for(int i=0; i<drugs.size(); i++)
		drug_cov[i]=from_map("drug_cov_" + to_string(i) + s, 0, 1);
	const double p=sum(drug_cov);
	if(p>1+1E-6)
		error_crit("drug coverages add up to more than 1, drug set " + to_string(n));
}

Human* Village::random_human(){
	Ran& ran=simulation->sim_ran;
	const int j=human_guide(ran);
	return village_het[j][ran_int(village_het[j].size(), ran)];
}
//***************************************************************************
//***************************************************************************

void Village::mass_birth(){
	const double t=simulation->event_manager.time_now();
	const int y=min(int(beta.size()-1), max(0, int(floor(t/dy + year0-init_demog_year))));

	Ran& ran=simulation->sim_ran;
	const int births=rpoisson(beta[y]*5.0/dy*total_pop(), ran);
	int b=0;
	while(b<births){
		Human* h=random_human();
		const double age=h->find_age(t)/dy;
		const int k=age<100.0 ? demog_age_guide(age/100.0) : alpha[y].size()-1;
		if(uniform(ran)<alpha[y][k]){
			h->die();
			b++;
		}
	}
}
void Village::update_mosq(const double dt){
	const double t=simulation->event_manager.time_now();
	num_updates++;
	const int nm=mosq_pop.size();
	for(int m=0; m<nm; m++)
		if(mosq_pop[m].M0>0)
			mosq_pop[m].update(dt);
}

void Village::update_foim(){
	if(det)
		parms.update_foim(false, parms.ode.time());
	else{
		const int nm=mosq_pop.size();
		const double t = simulation->event_manager.time_now();
		for(int j=0; j<village_het.size(); j++){
			const int vn=village_het[j].size();
			for(int i=0; i<vn; i++)
				village_het[j][i]->update_FOIM(false, t);
		}
		if (any_itn_irs) {
			
			for (int m = 0; m<nm; m++) {
				sum_repeat[m] = 0;
				sum_fail[m] = 0;
			}
			double sum_rel_biting_rate = 0;

			for (int j = 0; j<village_het.size(); j++) {
				for (int m = 0; m<nm; m++) {
					sum_repeat_j[m] = 0;
					sum_fail_j[m] = 0;
				}
				double sum_rel_biting_rate_j = 0;
				const int vn = village_het[j].size();
				for (int i = 0; i<vn; i++) {
					Human& h = *village_het[j][i];
					const double rel_biting_rate = simulation->human_parms.rel_biting_rate(h.find_age(t));
					for (int m = 0; m<nm; m++) {
						if (mosq_pop[m].M0>0) {
							sum_repeat_j[m] += rel_biting_rate*mosq_pop[m].prob_repeat(h.itn_irs, t);
							sum_fail_j[m] += rel_biting_rate*mosq_pop[m].prob_fail(h.itn_irs, t);
						}
					}
					sum_rel_biting_rate_j += rel_biting_rate;
				}
				const double het_exp_j = village_parms()->het_exp[j];
				for (int m = 0; m<nm; m++) {
					if (mosq_pop[m].M0>0) {
						sum_repeat[m] += sum_repeat_j[m]*het_exp_j;
						sum_fail[m] += sum_fail_j[m]*het_exp_j;
					}
				}
				sum_rel_biting_rate += sum_rel_biting_rate_j*het_exp_j;
			}

			for (int m = 0; m<nm; m++)
				if (mosq_pop[m].M0>0)
					mosq_pop[m].update_foim(sum_fail[m], sum_repeat[m], sum_rel_biting_rate);
		}
	}
}
int Village::total_pop() const{
	return N;
}
void Village::update_inc(const double age, const double prob_clin, const double prob_sev) {
	for (int l = 0; l<clin_inc.size(); l++)
		if (age>=clin_inc_age0[l] && age<clin_inc_age1[l])
			clin_inc[l] += prob_clin;
	for (int l = 0; l<inf_inc.size(); l++)
		if (age>=inf_inc_age0[l] && age<inf_inc_age1[l])
			inf_inc[l] += 1.0;
	for (int l = 0; l<sev_inc.size(); l++)
		if (age>=sev_inc_age0[l] && age<sev_inc_age1[l])
			sev_inc[l] += prob_sev;

}
void Village::output(const int index, const double time, const double output_interval){

	vector<vector<double> >& v=*simulation->to_output;
	vector<Output_var>& output_vars=simulation->output_vars;
	const int no=output_vars.size();
	const double t=simulation->event_manager.time_now(); 
	if(index>=0 && index<v[0].size()){
		if(det){
			vector<double> pred(parms.N_pop.size());
			parms.par_prev(parms.state, pred, false);
			vector<double> pred1(parms.N_pop.size());
			parms.par_prev(parms.state, pred1, true);
			v[0][index]=time;
			for(int l=0; l<no; l++){
				Output_var& o=output_vars[l];
				Output_var::Quantity q=o.quantity;
				if(q==Output_var::EIR)
					v[l+1][index]=parms.EIR[0]*dy;
				else if(q==Output_var::prev)
					v[l+1][index]=parms.mean_in_age_range(pred, o.age0, o.age1);
				else if(q==Output_var::pcr_prev)
					v[l+1][index]=parms.mean_in_age_range(pred1, o.age0, o.age1);
				else if(q==Output_var::clin_inc)
					v[l+1][index]=parms.mean_in_age_range(parms.clin_inc, o.age0, o.age1)*dy;
				else if(q==Output_var::inf_inc)
					v[l+1][index]=parms.mean_in_age_range(parms.inf_inc, o.age0, o.age1)*dy;
				else if(q==Output_var::sev_inc)
					v[l+1][index]=parms.mean_in_age_range(parms.sev_inc, o.age0, o.age1)*dy;
				else if(q==Output_var::prop)
					v[l+1][index]=parms.prop_in_age_range(o.age0, o.age1);
				else if(q==Output_var::itn_cov)
					v[l+1][index]=parms.itn_coverage;
			}
		}
		else{
			vector<double> N_age(no, 0);
			vector<double> prop(no, 0);
			const int n=village_het.size();
			for(int j=0; j<n; j++){
				Village_het& vh=village_het[j];
				const int nj=vh.size();
				for(int i=0; i<nj; i++){
					Human& h=*vh[i];
					const double age=h.find_age(t);

					for(int l=0; l<no; l++){
						Output_var& o=output_vars[l];
						Output_var::Quantity q=o.quantity;
						if(q != Output_var::EIR){
							if(age>=o.age0 && age<o.age1){
								N_age[l]++;
								if(q==Output_var::prev)
									prop[l] += h.prob_positive(h.find_age(t), false);
								else if(q==Output_var::pcr_prev)
									prop[l] += h.prob_positive(h.find_age(t), true);
								else if(q==Output_var::itn_cov && h.itn_irs.itn)
									prop[l]++;
							}
						}
					}
				}
			}

			v[0][index]=time;
			int inc_i=0;
			int inf_inc_i=0;
			int sev_inc_i=0;

			const double N1=total_pop();
			const double noy=dy/output_interval;
			for(int l=0; l<no; l++){
				Output_var::Quantity q=output_vars[l].quantity;
				
				if (q == Output_var::SM) {
					const int nm = mosq_pop.size();
					double SM = 0.0;
					for (int m = 0; m < nm; m++)
						if (mosq_pop[m].M0 > 0)
							SM += mosq_pop[m].S;
					v[l + 1][index] = SM;
				}

				if(q==Output_var::EIR)
					v[l+1][index]=(num_updates>0 ? EIR_count*dy/(village_parms()->update_interval*num_updates*N1) : -999);
				else if(q==Output_var::prev)
					v[l+1][index]=prop[l]/N_age[l];
				else if(q==Output_var::pcr_prev)
					v[l+1][index]=prop[l]/N_age[l];
				else if(q==Output_var::clin_inc)
					v[l+1][index]=clin_inc[inc_i++]*noy/N_age[l];
				else if(q==Output_var::inf_inc)
					v[l+1][index]=inf_inc[inf_inc_i++]*noy/N_age[l];
				else if(q==Output_var::sev_inc)
					v[l+1][index]=sev_inc[sev_inc_i++]*noy/N_age[l];
				else if(q==Output_var::prop)
					v[l+1][index]=N_age[l]/N1;
				else if(q==Output_var::itn_cov)
					v[l+1][index]=prop[l]/N_age[l];
			}
			int k=no+1;
			if(time>1E-5){
				v[k++][index]=num_itn;
				v[k++][index]=num_irs;
				v[k++][index]=num_trt;
				v[k++][index]=num_smc;
				v[k++][index]=num_mda_trt;
				v[k++][index]=num_mda_screen;
				v[k++][index]=num_vacc_doses;
			}
		}
	}

	EIR_count=0;
	num_updates=0;
	set_to_zero(clin_inc);
	set_to_zero(inf_inc);
	set_to_zero(sev_inc);
}

//***************************************************************************
//***************************************************************************

template<class F>
void Village::pulse(const double prob, const int target, const int j_repeated, F f, 
	const double age0, const double age1, const bool to_update_foim){
	/* target =
	0	use single-intervention correlations
	1	random
	2	most bitten : for 2 and 3, need large number (at least 10) of equally-sized heterogeneity groups
	3	least bitten
*/

	if(prob==0.0)
		return;
	const double t=simulation->event_manager.time_now();
	double mu0;
	if (target==0) {
		const double sigma = simulation->human_parms.sigma_repeated[j_repeated];
		mu0 = (prob<1.0 ? -sqrt(1+sigma*sigma)*invnorm(prob) : -1E12);
	}
	Ran& ran=simulation->sim_ran;
	for(int j=0; j<village_het.size(); j++){
		for(int i=0; i<village_het[j].size(); i++){
			Human* h=village_het[j][i];
			const double age=h->find_age(t);
			if(age<age0 || age>age1)
				continue;
			bool go=false;
			if(prob==1.0)
				go=true;
			else if(target==0)
				go=mu0 + h->muj_repeated[j_repeated] + rnorm(ran) < 0;
			else if(target==1)
				go=uniform(ran)< prob;
			else if(target==2){
				const int j1=village_het.size()-j-1;
				go=j1+1<prob*village_het.size() || (j1<prob*village_het.size() && uniform(ran)< prob*village_het.size()-j1);
			}
			else if(target==3)
				go=j+1<prob*village_het.size() || (j<prob*village_het.size() && uniform(ran)< prob*village_het.size()-j);
			else
				error_crit("target = " + to_string(target) +" not implemented");

			if(go)
				f(h);
		}
	}

	if(to_update_foim)
		update_foim();
}

//***************************************************************************
//***************************************************************************

Village_het::Village_het(){}

void Village_het::setup(Village* v, const int j){
	village=v;
	het_j=j;
}

//***************************************************************************
//***************************************************************************

double calc_immunity(double &I, double &I_time, double &I_boost_time, const double t, const double r,
					 const double u, const bool boost=true){
	const double dur=t-I_time;
	I_time=t;
	I *= exp(-dur*r);
	if (boost && (t-I_boost_time > u || I==0)) {
		I += 1.0;
		I_boost_time = t;
		return I-0.5;
	}
	else
		return I;
}

//***************************************************************************
//***************************************************************************

Human::Human(Village_het* vh) : village_het(vh), infection_state(S), 
		IC_time(-1), IB_time(-1), IV_time(-1), ID_time(-1),
		pevi(0), pev_efficacy(0.0), pev_ab0(0.0), pev_dose_times(0), num_vacc_doses(0),
		tbv(false), ipt(false), trt(-1), rel_c(0), trt_time(-1E20) {
	Village* village=village_het->village;
	hp = &village_het->village->simulation->human_parms;
	Ran& ran=village->simulation->sim_ran;

	move_event=Move_state_event(this);
	if(!village->demog){
		death=Death_event(this);
		death.schedule(&village->simulation->event_manager, -log(uniform())/hp->eta);
	}
	leave_itn_event=Leave_itn_event(this);
	epi_event=Epi_event(this);
	pev_dose_event=Pev_dose_event(this);

	FOIM_contribution.assign(village->mosq_pop.size(), 0);
	muj_repeated.resize(hp->sigma_repeated.size());
	rmnorm(hp->chol_repeated, muj_repeated, ran);
}
void Human::save_state(ofstream &out){
	Event_manager* em=&village_het->village->simulation->event_manager;
	save_event(out, em, move_event);
	save_event(out, em, death);
	save_event(out, em, leave_itn_event);
	save_event(out, em, epi_event);
	output_binary(epi_event.num_epi, out);
	::save_state(FOIM_contribution, out);
	output_binary(date_of_birth, out);
	output_binary(IC_M, out);
	output_binary(IC, out);
	output_binary(IC_time, out);
	output_binary(IC_boost_time, out);
	output_binary(IB, out);
	output_binary(IB_time, out);
	output_binary(IB_boost_time, out);
	output_binary(ID, out);
	output_binary(ID_time, out);
	output_binary(ID_boost_time, out);
	output_binary(IV_M, out);
	output_binary(IV, out);
	output_binary(IV_time, out);
	output_binary(IV_boost_time, out);
	itn_irs.save_state(out);
	output_binary(infection_state, out);
	output_binary(trt, out);
	output_binary(rel_c, out);               
	output_binary(trt_time, out);
}
void Human::use_state(ifstream &in){
	Event_manager* em=&village_het->village->simulation->event_manager;
	const double t0=time_now();
	use_event(in, em, move_event, t0);
	use_event(in, em, death, t0);
	use_event(in, em, leave_itn_event, t0);
	use_event(in, em, epi_event, t0);
	input_binary(epi_event.num_epi, in);
	::use_state(FOIM_contribution, in);
	input_binary(date_of_birth, in);
	input_binary(IC_M, in);
	input_binary(IC, in);
	input_binary(IC_time, in);
	input_binary(IC_boost_time, in);
	input_binary(IB, in);
	input_binary(IB_time, in);
	input_binary(IB_boost_time, in);
	input_binary(ID, in);
	input_binary(ID_time, in);
	input_binary(ID_boost_time, in);
	input_binary(IV_M, in);
	input_binary(IV, in);
	input_binary(IV_time, in);
	input_binary(IV_boost_time, in);
	itn_irs.use_state(in);
	input_binary(infection_state, in);
	input_binary(trt, in);
	input_binary(rel_c, in);               
	input_binary(trt_time, in);
}
void Human::save_event(ofstream &out, Event_manager* em, Event &ev){
	const bool scheduled=em->find(&ev);
	output_binary(scheduled, out);
	if(scheduled){
		const double t=ev.exe_time();
		output_binary(t, out);
	}
}
void Human::use_event(ifstream &in, Event_manager* em, Event &ev, const double t0){
	bool scheduled;
	input_binary(scheduled, in);
	if(scheduled){
		double t;
		input_binary(t, in);
		ev.schedule(em, t-t0);
	}
}
double Human::uniform() const{
	return ::uniform(village_het->village->simulation->sim_ran);
}
double Human::time_now() const{
	return village_het->village->simulation->event_manager.time_now();
}

void Human::set_state(Parms &parms, Human_state &state, const bool itn){
	Village& v=*village_het->village;
	Ran& ran=v.simulation->sim_ran;
	const double t=time_now();
	const int i=discrete_dev(parms.N_pop, ran);
	const double age=parms.age[i] + (i==parms.age.size()-1 ? 5.0*dy : parms.age[i+1]-parms.age[i])*uniform();

	date_of_birth=t-age;

	const int j=village_het->het_j;
	IC_M=parms.P_IC_M*parms.IC_20[(itn ? 1 : 0)][j];
	IV_M=parms.P_IV_M*parms.IV_20[(itn ? 1 : 0)][j];

	IC=state.IC_A[i];
	IB=state.IB[i];
	ID=state.ID[i];
	IV=state.IV[i];

	IC_time=t;
	IB_time=t;
	IC_boost_time=-1E20;
	IB_boost_time=-1E20;

	ID_time=t;
	IV_time=t;
	ID_boost_time=-1E20;
	IV_boost_time=-1E20;
	const int num_states=state.num_prev_states();
	for(int k=0; k<num_states; k++)
		hp->temp_prev[k]=state.state(k)[i];

	const double sump=sum(hp->temp_prev);
	for(int k=0; k<num_states; k++)
		hp->temp_prev[k] /= sump;

	const double u=uniform();
	double p=0;
	if(u< (p= p+hp->temp_prev[S]))
		infection_state=S;
	else if(u< (p= p+hp->temp_prev[T]))
		infection_state=T;
	else if(u< (p= p+hp->temp_prev[D]))
		infection_state=D;
	else if(u< (p= p+hp->temp_prev[A1]))
		infection_state=A1;
	else if(u< (p= p+hp->temp_prev[U]))
		infection_state=U;
	else
		infection_state=P;

	if(infection_state==T || infection_state==P){
		trt=discrete_dev(v.drug_cov, ran, false);
		const Drug& drug=v.drugs[trt];
		trt_time=t + drug.dur_T*log(uniform());
		if(infection_state==P){
			trt_time += drug.dur_P*log(uniform());
			infection_state=hp->next_state[T];
		}
	}

	handle_new_state(t);
}

void Human::handle_new_state(const double t){
	if(infection_state != S){
		Village& v = *village_het->village;

		if(infection_state==T)
			rel_c=v.drugs[trt].rc;

		double t_state;
		if(infection_state==P)
			t_state = max(v.drugs[trt].choose_dur(v.simulation->sim_ran) - (t - trt_time), 0.1);	
		else{
			double dur=1.0;
			if(infection_state==T)
				dur=v.drugs[trt].dur_T;		
			else if(infection_state==D)
				dur=hp->meanD;
			else if(infection_state==A1)
				dur=hp->meanA;
			else if(infection_state==U)
				dur=hp->meanU;
			t_state = -dur*log(uniform());
		}
		move_event.schedule(&v.simulation->event_manager, t_state);
	}
	update_FOIM(true, t);
}

void Human::update_FOIM(const bool new_state, const double t){
	if ((infection_state==S || infection_state==P) && !new_state)
		return;

	const double age = find_age(t);
	const double c_current = infectivity(age);
	if(c_current>0  || new_state){
		const int nm=FOIM_contribution.size();
		Village* village=village_het->village;
		copy(FOIM_contribution, village->old_FOIM_contribution);
		if(c_current>0.0){
			double F0 = c_current*hp->rel_biting_rate(age)*village->village_parms()->het_exp[village_het->het_j];
			if(tbv){
				const double d=hp->tbv_decay*(t-tbv_time);
				if(d>8)
					tbv=false;
				else
					F0 *= 1-hp->tbv_eff*exp(-d);
			}
			for(int m=0; m<nm; m++)
				FOIM_contribution[m]= (itn_irs.itn || itn_irs.irs) ? F0*(1-village->mosq_pop[m].prob_fail(itn_irs, t)) : F0;
		}
		else
			set_to_zero(FOIM_contribution);
		for(int m=0; m<nm; m++)
			village->mosq_pop[m].C += FOIM_contribution[m] - village->old_FOIM_contribution[m];
	}
}
void Human::move_state(){
	infection_state=hp->next_state[infection_state];
	handle_new_state(time_now());
}

void Human::possible_inf_bite(const double prob_bite, const double t){
	if(prob_bite>0.0){
		Village& v=*village_het->village;
		v.EIR_count += prob_bite; // EIR is for adults, so don't account for lower biting rate at younger age here
		double u=uniform();
		const double age=find_age(t);
		const double pba = prob_bite*hp->rel_biting_rate(age);

		if(u < pba){
			u /= pba;

			const double prob_inf =prob_infection(age);
			if(u < prob_inf){

				u /= prob_inf;
				const bool possible_clin=(infection_state!=T && infection_state!=D && infection_state!=P);
				const double prob_clin=prob_clinical_disease(age);
				const double prob_sev=prob_severe_disease(age)*(hp->separate_sev ? 1.0 : prob_clin);

				calc_immunity(ID, ID_time, ID_boost_time, t, hp->rd, hp->ud, true);

				if(possible_clin){
					v.update_inc(age, prob_clin, prob_sev);

					const Infection_state previous_state=infection_state;
					if(u < prob_clin){
						u /= prob_clin;
						const double f=sum(v.drug_cov);
						if(u<f){
							const int new_trt=discrete_dev(v.drug_cov, v.simulation->sim_ran, false);
							if(t>0.1)
								v.num_trt++;
							if(u<f*v.drugs[new_trt].efficacy){
								infection_state=T;
								trt=new_trt;
								trt_time=t;
							}
						}
						else
							infection_state=D;
					}
					else
						infection_state=A1;

					if(infection_state==A1 && previous_state==A1)
						return;
					if(previous_state != S)
						move_event.cancel(&v.simulation->event_manager);
					handle_new_state(t);
				}
			}
		}
	}
}

void Human::die(){

	Village& v=*village_het->village;

	Ran& ran=v.simulation->sim_ran;
	const double t=time_now();

	if(infection_state != S)
		move_event.cancel(&v.simulation->event_manager);

	// when someone dies, replace them in the same village with a new birth
	date_of_birth=t;
	if(!v.demog)
		death.schedule(&v.simulation->event_manager, -log(uniform())/hp->eta);

	for(int m=0; m<v.mosq_pop.size(); m++){
		v.mosq_pop[m].C -= FOIM_contribution[m];
		FOIM_contribution[m]=0;
	}

	infection_state=S;
	trt_time=-1E20;
	trt=-1;

	bool suitable=false;
	Human* h=0;
	int num_tries=0;
	while(!suitable){
		h= (num_tries++>100) ? v.random_human() : (*village_het)[ran_int(village_het->size(), ran)];
		const double a=h->find_age(t);
		suitable=(num_tries<200 ? (a>15*dy && a<35*dy) : true);
	}

	IC_M=hp->P_IC_M*h->IC;
	IV_M=hp->P_IV_M*h->IV;
	IC=0;
	IC_time=t;
	IC_boost_time=-1E20;
	IB=0;
	IB_time=t;
	IB_boost_time=-1E20;

	IV=0;
	IV_time=t;
	IV_boost_time=-1E20;
	ID=0;
	ID_time=t;
	ID_boost_time=-1E20;

	if(!hp->retain_itn_irs){
		if(hp->leave_itn_on && itn_irs.itn)
			leave_itn_event.cancel(&v.simulation->event_manager);
		itn_irs=ITN_IRS();
	}
	if(v.prob_itn_birth>0){
		const double u=uniform();
		if(u<v.prob_itn_birth){
			give_itn();
			if(u<v.prob_itn_birth_adult)
				h->give_itn();
		}
	}

	pevi=0;
	pev_dose_times.clear();
	pev_efficacy=0;
	pev_ab0=0;

	pev_dose_event.cancel(&v.simulation->event_manager);

	tbv=false;
	ipt=false;
	num_vacc_doses=0;
	
	epi_event.cancel(&v.simulation->event_manager);
	epi_event.num_epi=0;
	epi_event.schedule(&v.simulation->event_manager, v.simulation->epi_ages[0] + 10E-6*uniform());

}

double Human::infectivity(const double age){
	if(infection_state==S || infection_state==P)
		return 0.0;
	else
		return hp->infectivity(infection_state, prob_detection(age), rel_c);
}

double Human::prob_positive(const double age, const bool pcr){
	if(infection_state==S || infection_state==P)
		return 0.0;
	else if(infection_state==D || infection_state==T)
		return 1.0;
	else
		return hp->prob_positive(infection_state, prob_detection(age), pcr);
}

//***************************************************************************

void Human::give_itn(){
	Village& v=*village_het->village;
	v.num_itn++;
	v.any_itn_irs=true;
	Event_manager* event_manager=&v.simulation->event_manager;
	if(hp->leave_itn_on){
		if(itn_irs.itn)
			leave_itn_event.cancel(event_manager);
		leave_itn_event.schedule(event_manager, -hp->leave_itn_time*log(uniform()));
	}
	itn_irs.itn=true;
	itn_irs.itn_time=time_now();
}

void Human::leave_itn(){
	itn_irs.itn=false;
	update_FOIM(false, time_now());
}

void Human::give_irs(){
	Village& v=*village_het->village;
	v.num_irs++;
	v.any_itn_irs=true;
	itn_irs.irs=true;
	itn_irs.irs_time=time_now();
}

void Human::pev_dose() {
	Village& v = *village_het->village;
	num_vacc_doses++;
	v.num_vacc_doses++;
	if (num_vacc_doses<pev_dose_times.size())
		pev_dose_event.schedule(&v.simulation->event_manager,
			pev_dose_times[num_vacc_doses]-pev_dose_times[num_vacc_doses-1]);

	if (num_vacc_doses==1 || num_vacc_doses==4) {
		Ran& ran = v.simulation->sim_ran;
		Pev& pev = hp->pevs[pevi-1];
		const bool boost = num_vacc_doses==4;
		if (pev.ab_model) {
			const double mu = boost ? pev.ab0_boost_mu : pev.ab0_mu;
			const double sigma = boost ? pev.ab0_boost_sigma : pev.ab0_sigma;
			const double rho_mu = boost ? pev.rho_boost_mu : pev.rho_mu;
			const double rho_sigma = boost ? pev.rho_boost_sigma : pev.rho_sigma;
			pev_ab0 = exp(mu + sigma*rnorm(ran));
			pev_d1 = exp(pev.d1_mu + pev.d1_sigma*rnorm(ran));
			pev_d2 = exp(pev.d2_mu + pev.d2_sigma*rnorm(ran));
			pev_rho = invlogit(rho_mu + rho_sigma*rnorm(ran));
		}
		else {
			const double e = boost ? pev.boost_efficacy : pev.efficacy;
			const double pev_c = pev.leakiness;
			if (pev_c<=0.1) // all or nothing
				pev_efficacy = uniform()<e ? 1 : 0;
			else if (pev_c>0.1 && pev_c<10000)
				pev_efficacy = rbeta(pev_c*e, pev_c*(1-e), ran);
			else if (pev_c>=10000)	 // wholly leaky
				pev_efficacy = e;
		}
	}
}

void Human::give_pev(const int vaccine, const double boost_cov, const vector<double> &dose_intervals){
	pevi = vaccine;
	num_vacc_doses = 0;
	const double t = time_now();
	pev_dose_times.resize(dose_intervals.size()>3 && uniform()<boost_cov ? dose_intervals.size() : 3);
	for (int i = 0; i<pev_dose_times.size(); i++)
		pev_dose_times[i] = t + dose_intervals[i]-dose_intervals[0];
	
	pev_dose_event.cancel(&village_het->village->simulation->event_manager);
	pev_dose();
}
void Human::give_tbv(){
	tbv=true;
	tbv_time=time_now();
}

void Human::give_epi(const int num_epi){
	Village& v=*village_het->village;

	if(num_epi==0 && (v.prob_pev_epi>0 || v.prob_ipt_epi>0)){
		Ran& ran=v.simulation->sim_ran;
		const int j_repeated=4;
		bool go[2];
		double prob[2]={v.prob_pev_epi, v.prob_ipt_epi};
		const double z=rnorm(ran);
		for(int k=0; k<2; k++){
			if (prob[k]==0.0)
				go[k] = false;
			else if (prob[k]==1.0)
				go[k] = true;
			else {
				const double sigma = hp->sigma_repeated[j_repeated];
				const double mu0 = -sqrt(1+sigma*sigma)*invnorm(prob[k]);
				go[k] = mu0 + muj_repeated[j_repeated] + z < 0;
			}
		}
		if (go[0])
			give_pev(1, hp->pev_epi_boost_coverage, v.simulation->epi_ages);
		ipt=go[1];
	}
	

	if(ipt)
		treat(v.ipt_drug);
}
void Human::treat(const int new_trt){
	// this is only for presumptive treatment - MDA, MSAT, IPT etc.
	Village& v=*village_het->village;
	const Drug& drug=v.drugs[new_trt];
	if(uniform() < drug.efficacy){
		if(infection_state != S)
			move_event.cancel(&v.simulation->event_manager);

		const double t = time_now();
		trt_time=t;
		trt=new_trt;

		infection_state = (infection_state==D || (infection_state==A1 && uniform()<prob_detection(find_age(t))))
			? T : (drug.p_protect ? S : P);
		handle_new_state(t);
	}
}

//***************************************************************************  
//***************************************************************************

double Human::find_age() const{
	return find_age(time_now());
}
double Human::find_age(const double t) const {
	return t - date_of_birth;
}
int Human::find_age_group(const double age) const{
	return village_het->village->parms.find_age_group(age);
}
//***************************************************************************

double Human::prob_infection(const double age){
	const double t=time_now();
	const double pp=(infection_state==S && trt>=0 
		? village_het->village->drugs[trt].prob_infected(t-trt_time, age) 
		: (infection_state==P  ? 0.0 : 1.0));
	if(pp>0.999 && infection_state==S)
		trt=-1;

	const double IB1=calc_immunity(IB, IB_time, IB_boost_time, t, hp->rb, hp->ub, true);
	if (pp<1E-6)
		return 0.0;

	double p=pp*hp->prob_infection(IB1);
	if(pevi>0)
		p *= 1.0 - hp->pevs[pevi-1].find_efficacy(*this, t-pev_dose_times[num_vacc_doses-1]);
	return p;
}

//***************************************************************************

double Human::prob_clinical_disease(const double age){
	const double t=time_now();
	const double IC1=calc_immunity(IC, IC_time, IC_boost_time, t, hp->rc, hp->uc, true);
	return hp->prob_clinical_disease(IC1 + IC_M*exp(-age/hp->dm));
}
double Human::prob_severe_disease(const double age){
	const double t=time_now();
	const double IV1=calc_immunity(IV, IV_time, IV_boost_time, t, hp->rv, hp->uv, true);
	return hp->prob_severe_disease(IV1 + (age>8*hp->dvm ? 0 : IV_M*exp(-age/hp->dvm)), find_age_group(age));
}
double Human::prob_detection(const double age){
	const double t=time_now();
	const double ID1=calc_immunity(ID, ID_time, ID_boost_time, t, hp->rd, hp->ud, false);
	return hp->prob_detection(ID1, find_age_group(age));
}
//***************************************************************************

double Human_parms::rel_biting_rate(const double age){
	return rel_biting_rate_table[min(int(age*br0), Tn-1)];
}
double Human_parms::prob_infection(const double IB1){
	return prob_inf_table[min(int(IB1*IB00), Tn-1)];
}
double Human_parms::prob_clinical_disease(const double IC1){
	return prob_clin_table[min(int(IC1*IC00), Tn-1)];
}
double Human_parms::prob_detection(const double I, const int ai){
	return (1-dmin)/(1+prob_det_table[min(int(I*ID00), Tn-1)]*fD[ai]) + dmin;
}
double Human_parms::prob_severe_disease(const double I, const int ai){
	return theta0*((1-theta1)/(1+prob_sev_table[min(int(I*IV00), Tn-1)]*fV[ai]) + theta1);
}
double Human_parms::infectivity(const Human::Infection_state infection_state, const double prob_det, const double rel_c){
	if(infection_state==Human::D)
		return cD;
	else if(infection_state==Human::T)
		return cD*rel_c;
	else if(infection_state==Human::U)
		return cU;
	else if(infection_state==Human::A1)
		return cA_table[min(int(prob_det*Tn), Tn-1)];
	else
		return 0.0;
}
double Human_parms::find_pcr_par(const double q){
	return pcrA_table[min(int(q*Tn), Tn-1)];
}

double Human_parms::prob_positive(const Human::Infection_state infection_state, const double prob_det, const bool pcr){
	if(infection_state==Human::D || infection_state==Human::T)
		return 1.0;
	else if(infection_state==Human::U)
		return pcr ? pcrU_table[min(int(prob_det*Tn), Tn-1)] : 0.0;
	else if(infection_state==Human::A1)
		return pcr ? pcrA_table[min(int(prob_det*Tn), Tn-1)] : prob_det;
	else
		return 0.0;

}

//***************************************************************************  
//***************************************************************************

