//***************************************************************************
//***************************************************************************

Mosq_pop::Mosq_pop(Village* v) : village(v), parms(0), first(true), w_bar(1), z_bar(0), C(0), 
	tau1(0.69), rel_M(1), seasonal_table(365), birth_index(0), larvicide(false){ }

void Mosq_pop::setup(Parms* parms1, const string &name, const int N0, const string &pop_file){
	
	ran.setup(village->simulation->sim_ran.ran64());
	species = name;

	N=N0;
	parms=parms1;

	year0=from_map_int("year0", 1900, 2100, 2014);
	if(pop_r.size()==0){
		growing_pop= village->demog ? from_map_bool("growing_pop", false) : false;
		if(growing_pop){
			pop_r.clear();
			vector<int> pop_years;
			import_file(pop_file, pop_years, pop_r);
			init_pop_year=pop_years.front();
		}
		else{
			pop_r.assign(1, 1.0);
			init_pop_year=2010;
		}
	}

	f		=	from_map(species+"_biting_rate", 0, 1/tau1);
	mu0		=	from_map(species+"_mu0");
	Q0		=	from_map(species+"_Q0", 0, 1);
	Q_in	=	from_map(species+"_Q_in", 0, 1);
	Q_bed	=	from_map(species+"_Q_bed", 0, 1);
		
	if(in_map("itn_dnw_"+species) && !in_map("itn_kill_"+species))
		itn_irs_2016=true;
	else if(!in_map("itn_dnw_"+species) && in_map("itn_kill_"+species))
		itn_irs_2016=false;
	else if(!in_map("itn_dnw_"+species) && !in_map("itn_kill_"+species))
		error_crit("Need to specify ITN/IRS parameters, either pre- or post-2016 set");
	else
		itn_irs_2016=from_map_bool("itn_irs_2016", false);

	irs_ellie=from_map_bool("irs_ellie", false);

	if(irs_ellie && itn_irs_2016)
		error_crit("Can't have irs_ellie=1 and itn_irs_2016=1");
	
	if(itn_irs_2016){
		itn_rn = in_map("itn_rn_"+species) ? from_map("itn_rn_"+species,0,1) : from_map("itn_re_"+species,0,1);
		itn_rnw = in_map("itn_rnw_"+species) ? from_map("itn_rnw_"+species,0,1) : from_map("itn_rp_"+species,0,1);
		itn_dnw = from_map("itn_dnw_"+species,0,1);
		itn_dnf = from_map("itn_dnf_"+species,0,1,0);
		irs_ri = in_map("irs_ri_"+species) ? from_map("irs_ri_"+species,0,1) : from_map("irs_repel_"+species, 0, 1);
		irs_diw = from_map("irs_diw_"+species,0,1,0);
		irs_dif =in_map("irs_dif_"+species) ?  from_map("irs_dif_"+species,0,1) : from_map(species+"_endophily", 0, 1); 
		irs_riw = from_map("irs_riw_"+species,0,1,0);
	}
	else{
		itn_kill	=	from_map("itn_kill_"+species, 0, 1);
		itn_repel	=	from_map("itn_repel_"+species, 0, 1);
		itn_repel_min=	from_map("itn_repel_min_"+species, 0, 1);
		
		if(irs_ellie){
			irs_decay_mort1 = from_map("irs_decay_mort1", -1E20);
			irs_decay_mort2 = from_map("irs_decay_mort2", -1E20);
			irs_decay_succ1 = from_map("irs_decay_succ1", -1E20);
			irs_decay_succ2 = from_map("irs_decay_succ2", -1E20);
			irs_decay_det1 = from_map("irs_decay_det1", -1E20);
			irs_decay_det2 = from_map("irs_decay_det2", -1E20);
			irs_k0 = from_map("irs_k0", 0, 1);
			irs_succ_k0 = from_map_bool("irs_succ_k0", true);
		}
		else{
			irs_kill	=	from_map("irs_kill_"+species, 0, 1);
			irs_repel	=	from_map("irs_repel_"+species, 0, 1);
			endophily	=	from_map(species+"_endophily", 0, 1);
		}
	}
	
	feedback=from_map_bool("larval_feedback");
	stochastic=from_map_bool("stochastic_mosq");
	itn_decay_rate=log(2.0)/(dy*from_map("itn_half_life"));
	irs_decay_rate=log(2.0)/(dy*from_map("irs_half_life"));
	itn_decay_shape=in_map("itn_decay_shape") ? from_map("itn_decay_shape") : 1.0;
	irs_decay_shape=in_map("irs_decay_shape") ? from_map("irs_decay_shape") : 1.0;

	M0=from_map("prop_"+species)*from_map("total_M")*N*from_map("rel_mv");

	bool use_fourier;
	if(in_map("seasonal_a0")){
		if(in_map(species+"_k"))
			error_crit("Can't have both Fourier and non-Fourier seasonal parameters in input file");
		use_fourier=true;
	}
	else
		use_fourier=false;

	const double two_pi = 6.283185307179586;
	if(use_fourier){
		int num_fourier=0;
		while(in_map("seasonal_a"+to_string(num_fourier+1)))
			num_fourier++;
		const double a0=from_map("seasonal_a0");
		vector<double> a(num_fourier);
		vector<double> b(num_fourier);
		for(int j=0; j<num_fourier; j++){
			a[j]=from_map("seasonal_a"+to_string(j+1), -1E20);
			b[j]=from_map("seasonal_b"+to_string(j+1), -1E20);
		}
		double seasonal_mean = 0.0;
		for (int i = 0; i<365; i++) {
			const double x = two_pi*(i+0.5)/dy;
			double y = a0;
			for (int j = 0; j<num_fourier; j++)
				y += a[j]*cos(x*(j+1)) + b[j]*sin(x*(j+1));
			seasonal_table[i] = y;
			seasonal_mean += seasonal_table[i];
		}
		seasonal_mean /= 365.0;
		for (int i = 0; i<365; i++)
			seasonal_table[i] = max(seasonal_table[i]/seasonal_mean, 0.01);
	}
	else{
		const double k=from_map(species+"_k");
		const double u=from_map(species+"_u", 0, 1);
		const double c=from_map(species+"_c", 0, 1);
		double seasonal_mean = 0.0;
		for (int i = 0; i<365; i++) {
			const double x = two_pi*((i+0.5)/dy - u);
			seasonal_table[i] = pow(0.5*(1+cos(x)), k);
			seasonal_mean += seasonal_table[i];
		}
		seasonal_mean /= 365.0;
		for (int i = 0; i<365; i++)
			seasonal_table[i] = c + (1-c)*seasonal_table[i]/seasonal_mean;

	}


	if(feedback){
		mue=from_map("larval_mue");
		mul=from_map("larval_mul");

		de=from_map("larval_de");
		dl=from_map("larval_dl");
		dp=from_map("larval_dp");
		const double beta0=from_map("larval_beta");
		const double x=exp(-mu0/f);
		eov=beta0*(1-x)/(x*mu0);
		mup=from_map("larval_mup");
		gamma=from_map("larval_gamma");
		b_lambda=(gamma*mul/mue-de/dl+(gamma-1)*mul*de);
	}

	latmosq=from_map("latmosq");
	latgam=parms->latgam;

	omega=parms->omega;
	D=omega*N;

	tau2=1/f-tau1;
	av=f*Q0;

	mu=mu0;
	p1=exp(-mu0*tau1);
	p2=exp(-mu0*tau2);

	Q=Q0;
	PM=exp(-mu*latmosq);


	find_init();

	infv=0;
	lag_infv.setlag(latmosq);
}

void Mosq_pop::save_state(ofstream &out){
	output_binary(infv, out);
	lag_infv.save_state(out);
	ran.save_state(out);
	output_binary(first, out);
	output_binary(FOIM, out);
	output_binary(FOIM_end, out);
	lag_FOIM.save_state(out);
	lag_eir.save_state(out);		
	lag_inc_E.save_state(out);
	output_binary(FOIM_init, out);
	output_binary(C, out);
	output_binary(rel_M, out);
	output_binary(K0, out);
	output_binary(M0, out);
	output_binary(S, out);
	output_binary(E, out);
	output_binary(I, out);
	output_binary(EL, out);
	output_binary(LL, out);
	output_binary(PL, out);
	output_binary(p1, out);
	output_binary(p2, out);
	output_binary(av, out);
	output_binary(mu, out);
	output_binary(PM, out);
	output_binary(f, out);
	output_binary(Q, out);
	output_binary(D, out);
	output_binary(w_bar, out);
	output_binary(z_bar, out);
}

void Mosq_pop::use_state(ifstream &in){
	input_binary(infv, in);
	lag_infv.use_state(in);
	ran.use_state(in);
	input_binary(first, in);
	input_binary(FOIM, in);
	input_binary(FOIM_end, in);
	lag_FOIM.use_state(in);
	lag_eir.use_state(in);		
	lag_inc_E.use_state(in);
	input_binary(FOIM_init, in);
	input_binary(C, in);
	input_binary(rel_M, in);
	input_binary(K0, in);
	input_binary(M0, in);
	input_binary(S, in);
	input_binary(E, in);
	input_binary(I, in);
	input_binary(EL, in);
	input_binary(LL, in);
	input_binary(PL, in);
	input_binary(p1, in);
	input_binary(p2, in);
	input_binary(av, in);
	input_binary(mu, in);
	input_binary(PM, in);
	input_binary(f, in);
	input_binary(Q, in);
	input_binary(D, in);
	input_binary(w_bar, in);
	input_binary(z_bar, in);
}

void Mosq_pop::reset(){

	f=1/(tau1+tau2);
	mu=mu0;
	D=omega*N;
	av=f*Q0;
	p1=exp(-mu0*tau1);
	p2=exp(-mu0*tau2);
	Q=Q0;
	PM=exp(-mu*latmosq);

	w_bar=1;
	z_bar=0;
	C=0;
	find_init();
}

void Mosq_pop::change_itn(const int n){
	const string suffix="_" + to_string(n);
	if(itn_irs_2016){
		itn_rn=from_map_suffix("itn_rn_"+species,suffix,0,1);
		itn_rnw=from_map_suffix("itn_rnw_"+species,suffix,0,1);
		itn_dnw=from_map_suffix("itn_dnw_"+species,suffix,0,1);
		itn_dnf=from_map_suffix("itn_dnf_"+species,suffix,0,1, 0.0);
	}
	else{
		itn_kill=from_map_suffix("itn_kill_"+species, suffix, 0, 1);
		itn_repel=from_map_suffix("itn_repel_"+species, suffix, 0, 1);
		itn_repel_min=from_map_suffix("itn_repel_min_"+species, suffix, 0, 1);
	}
	itn_decay_rate=log(2.0)/(dy*from_map_suffix("itn_half_life", suffix));
}

void Mosq_pop::change_irs(const int n){
	const string suffix="_" + to_string(n);
	if(itn_irs_2016){
		irs_ri = from_map_suffix("irs_ri_"+species, suffix, 0,1);
		irs_diw = from_map_suffix("irs_diw_"+species, suffix,0, 1, 0);
		irs_dif =from_map_suffix("irs_dif_"+species, suffix,0,1); 
		irs_riw = from_map_suffix("irs_riw_"+species, suffix,0,1,0);
	}
	else{
		if(irs_ellie){
			irs_decay_mort1 = from_map_suffix("irs_decay_mort1", suffix, -1E20);
			irs_decay_mort2 = from_map_suffix("irs_decay_mort2", suffix, -1E20);
			irs_decay_succ1 = from_map_suffix("irs_decay_succ1", suffix, -1E20);
			irs_decay_succ2 = from_map_suffix("irs_decay_succ2", suffix, -1E20);
			irs_decay_det1 = from_map_suffix("irs_decay_det1", suffix, -1E20);
			irs_decay_det2 = from_map_suffix("irs_decay_det2", suffix, -1E20);
			irs_k0 = from_map_suffix("irs_k0", suffix, 0, 1);
			irs_succ_k0 = from_map_bool_suffix("irs_succ_k0", suffix);
		}
		else{
			irs_kill	=	from_map_suffix("irs_kill_"+species, suffix, 0, 1);
			irs_repel	=	from_map_suffix("irs_repel_"+species, suffix, 0, 1);
		}
	}
	irs_decay_rate=log(2.0)/(dy*from_map_suffix("irs_half_life", suffix));
}

void Mosq_pop::change_larvicide(const double p, const int num) {
	if (p==0.0)
		larvicide = false;
	else {
		larvicide = true;
		larval_time = village->simulation->event_manager.time_now();
		const string suffix = num==1 ? "_" + species : "_" + species + "_" + to_string(num);
		deprec_param = from_map_suffix("deprec_param", suffix, -1E20, 1E20);
		larvi_min = from_map_suffix("larvi_min", suffix, 0.0, 1.0);
	}
}

//***************************************************************************
// used only in deterministic model
void Mosq_pop::update_history(const double t){
	if(parms->run_type<2)
		lag_infv.update(t, infv);
	else{
		if(b_rate_t.size()==0 || t>b_rate_t.back()){
			b_rate.push_back(0.5*PL/dp);
			b_rate_t.push_back(t);
		}
	}
}
void Mosq_pop::update_det(const Mosq_state &y,  Mosq_state &dydt, const double t){
	/* run_type
	0	run full model (with larval stages if feedback == 1)
	1	run transmission model (human states and adult mosquito states), using stored mosquito birth rate
	2	run just larval and adult mosquito model and store adult mosquito birth rate
	*/
	double births;
	if(parms->run_type!=1)
		births=(feedback ? 0.5*y.PL/dp : rel_M*mu0*M0*seasonal_curve(t));
	else{
		double t1=t;
		while(t1<b_rate_t[0])
			t1 += dy;
		while(birth_index>0 && b_rate_t[birth_index]>t1)
			birth_index--;
		while(birth_index<b_rate_t.size()-1 && b_rate_t[birth_index+1]<t1)
			birth_index++;
		if(birth_index==b_rate_t.size()-1)
			births=b_rate[birth_index];
		else
			births=b_rate[birth_index] + (b_rate[birth_index+1]-b_rate[birth_index])*(t1-b_rate_t[birth_index])/(b_rate_t[birth_index+1]-b_rate_t[birth_index]);
	}

	if(parms->run_type==2)
		dydt.S +=  - mu*y.S + births;
	else{
		const double FOIM1=Q*f*C/D;
		infv=max(0.0, FOIM1*y.S);
		const double incv=max(0.0, lag_infv(t)*PM);
		dydt.S += -infv - mu*y.S + births;
		dydt.E += infv-incv - mu*y.E;
		dydt.I += incv - mu*y.I;
	}
	I=y.I;
	S=y.S;
	E=y.E;
	if(feedback && parms->run_type!=1){
		const double M=S + (parms->run_type==2 ? 0 : I + E); 
		const double K=K0*seasonal_curve(t);
		find_beta();
		const double max_mu=5;
		dydt.EL += beta*M - min(max_mu, mue*(1+(y.EL+y.LL)/K))*y.EL-y.EL/de;
		dydt.LL += y.EL/de - min(max_mu, mul*(1+gamma*(y.EL+y.LL)/K))*y.LL - y.LL/dl;
		dydt.PL += y.LL/dl - mup*y.PL - y.PL/dp;
		LL=y.LL;
		EL=y.EL;
		PL=y.PL;
	}
}
void Mosq_pop::find_equilibrium_det(Mosq_state &y, const double mean_c){
	const double foiv=mean_c*av;
	double Sp, Ep, Ip;
	if(foiv<1E-11 || parms->run_type==2){
		Sp=1;
		Ep=0;
		Ip=0;
	}
	else{
		Ip=foiv*PM/(foiv+mu);
		Sp=mu*Ip/(foiv*PM);
		Ep=1-Sp-Ip;
	}
	y.S=Sp*M0*rel_M;
	y.E=Ep*M0*rel_M;
	y.I=Ip*M0*rel_M;
	infv=foiv*y.S;
	S=y.S;
	I=y.I;
	E=y.E;

	if(feedback){
		find_init();
		y.EL=EL;
		y.LL=LL;
		y.PL=PL;
	}
}
//***************************************************************************

void Mosq_pop::find_beta(){
	const double x=exp(-mu/f);
	beta=eov*mu*x/(1-x);
}

void Mosq_pop::find_init(){
	const double M=rel_M*M0;
	if(M>0 && feedback){
		find_beta();
		const double lambda=-0.5*b_lambda + sqrt(0.25*b_lambda*b_lambda + gamma*beta*mul*de/(2*mue*mu*dl*(1+dp*mup)));
		K0=2*M*dl*mu*(1+dp*mup)*gamma*(lambda+1)/(lambda/(mul*de)-1/(mul*dl)-1);
		PL=2*dp*mu*M;
		LL=dl*(mup+1/dp)*PL;
		EL=(LL/dl + mul*LL*(1+gamma*LL/K0))/(1/de-mul*gamma*LL/K0);
	}
	else{
		K0=0;
		PL=0;
		LL=0;
		EL=0;
	}
}

void Mosq_pop::set_init_indiv(Mosq_state &y, const double t){

	FOIM=Q*f*C/D;
	FOIM_end=FOIM;
	FOIM_init=FOIM;

	I=y.I;
	S=y.S;
	E=y.E;
	if(feedback){
		EL=y.EL;
		LL=y.LL;
		PL=y.PL;
	}

	if(stochastic){
		I=floor(I+0.5);
		E=floor(E+0.5);
		S=floor(S+0.5);
	}

	const double update_step=1.0/village->simulation->mosq_update_dy;
	const double inc_E=floor(S*FOIM*update_step+0.5);
	lag_FOIM.setup(latgam, update_step, FOIM);
	lag_eir.setup(parms->dur_E, update_step, find_eir(t));		
	lag_inc_E.setup(latmosq, update_step, inc_E);
}

double Mosq_pop::itn_decay(const double t_ago){
	const double log2=0.6931471805599453;
	return itn_decay_shape==1.0 ? exp(-itn_decay_rate*t_ago) : 
			exp(-log2*pow(itn_decay_rate*t_ago/log2, itn_decay_shape));
}
double Mosq_pop::irs_decay(const double t_ago){
	const double log2=0.6931471805599453;
	return irs_decay_shape==1.0 ? exp(-irs_decay_rate*t_ago) : 
			exp(-log2*pow(irs_decay_rate*t_ago/log2, irs_decay_shape));
}

/*
For each of the following three functions, the probability (given that the attempt is on a human) is calculated as:
	(1-Q_in)	x	Pr(outcome if feeding attempt is outdoors) 
 +	Q_bed		x	Pr(outcome with ITN and IRS protection at the level this person has) 
 + (Q_in-Q_bed) x	Pr(outcome with just IRS protection at the level this person has)
 */

double Mosq_pop::prob_fail(const ITN_IRS &itn_irs, const double t){
	if(!itn_irs.itn && !itn_irs.irs)
		return 0.0;

	const double r_itn=itn_irs.itn ? itn_decay(t-itn_irs.itn_time) : 0.0;
	const double r_irs=itn_irs.irs ? irs_decay(t-itn_irs.irs_time) : 0.0;
	double w_irs=1.0;
	double w_itn_irs=1.0;
	if(itn_irs_2016){
		const double rn=itn_rn*r_itn;
		const double rnw=itn_rnw*r_itn;
		const double dnw=itn_dnw*r_itn;
		const double dnf=itn_dnf*r_itn;
		const double ri=irs_ri*r_irs;
		const double diw=irs_diw*r_irs;
		const double dif=irs_dif*r_irs;
		const double riw=irs_riw*r_irs;
		w_irs=(1-ri)*(1-diw-riw)*(1-dif);
		w_itn_irs=w_irs*(1-rn)*(1-rnw-dnw)*(1-dnf);
	}
	else{
		const double r1=(itn_irs.itn ? itn_repel_min + (itn_repel-itn_repel_min)*r_itn : 0);
		const double d1=itn_kill*r_itn;
		double r2, d2;
		if(irs_ellie && itn_irs.irs){
			const double ti=t-itn_irs.irs_time;
			const double det_hut_irs = 1 / (1 + exp(-(irs_decay_det1 + irs_decay_det2*ti)));
			const double succ_hut_irs = (irs_succ_k0 ? irs_k0 : 1.0) / (1 + exp(-(irs_decay_succ1 + irs_decay_succ2*ti)));
			const double mort_hut_irs = 1 / (1 + exp(-(irs_decay_mort1 + irs_decay_mort2*ti)));
			const double rep_hut_irs = 1 -  succ_hut_irs -  mort_hut_irs;
			const double kp_irs = (1 - det_hut_irs)*succ_hut_irs;
			const double jp_irs = (1 - det_hut_irs)*rep_hut_irs + det_hut_irs;
			const double lp_irs = (1 - det_hut_irs)*mort_hut_irs;
			r2=(1 - kp_irs/irs_k0)*(jp_irs/(lp_irs + jp_irs));
			d2=(1 - kp_irs/irs_k0)*(lp_irs/(lp_irs + jp_irs))/(1-r2);
		}
		else if(itn_irs.irs){
			r2=irs_repel*r_irs;
			d2=endophily*r_irs;
		}
		else{
			r2=0.0;
			d2=0.0;
		}
		w_irs=(1-r2)*(1-d2);
		w_itn_irs=w_irs*(1-d1-r1);
	}
	return Q_bed*(1-w_itn_irs) + (Q_in-Q_bed)*(1-w_irs);
}
double Mosq_pop::prob_repeat(const ITN_IRS &itn_irs, const double t){
	if(!itn_irs.itn && !itn_irs.irs)
		return 0.0;
	const double r_itn=itn_irs.itn ? itn_decay(t-itn_irs.itn_time) : 0.0;
	const double r_irs=itn_irs.irs ? irs_decay(t-itn_irs.irs_time) : 0.0;

	double z_irs=0;
	double z_itn_irs=0;
	if(itn_irs_2016){
		const double rn=itn_rn*r_itn;
		const double rnw=itn_rnw*r_itn;
		const double ri=irs_ri*r_irs;
		const double diw=irs_diw*r_irs;
		const double riw=irs_riw*r_irs;
		z_irs=ri + (1-ri)*riw;
		z_itn_irs=ri + (1-ri)*(rn + (1-rn)*(riw + (1-diw-riw)*rnw));
	}
	else{
		const double r1=(itn_irs.itn ? itn_repel_min + (itn_repel-itn_repel_min)*r_itn : 0);
		double r2;
		if(irs_ellie && itn_irs.irs){
			const double ti=t-itn_irs.irs_time;
			const double det_hut_irs = 1 / (1 + exp(-(irs_decay_det1 + irs_decay_det2*ti)));
			const double succ_hut_irs = (irs_succ_k0 ? irs_k0 : 1.0) / (1 + exp(-(irs_decay_succ1 + irs_decay_succ2*ti)));
			const double mort_hut_irs = 1 / (1 + exp(-(irs_decay_mort1 + irs_decay_mort2*ti)));
			const double rep_hut_irs = 1 -  succ_hut_irs -  mort_hut_irs;
			const double kp_irs = (1 - det_hut_irs)*succ_hut_irs;
			const double jp_irs = (1 - det_hut_irs)*rep_hut_irs + det_hut_irs;
			const double lp_irs = (1 - det_hut_irs)*mort_hut_irs;
			r2=(1 - kp_irs/irs_k0)*(jp_irs/(lp_irs + jp_irs));
		}
		else if(itn_irs.irs){
			r2=irs_repel*r_irs;
		}
		else
			r2=0.0;
		z_irs=r2;
		z_itn_irs=r2+(1-r2)*r1;
	}
	return Q_bed*z_itn_irs + (Q_in-Q_bed)*z_irs;
} 
double Mosq_pop::prob_bite(const ITN_IRS &itn_irs, const double t){
	if(!itn_irs.itn && !itn_irs.irs)
		return 1.0;

	const double r_itn=itn_irs.itn ? itn_decay(t-itn_irs.itn_time) : 0.0;
	const double r_irs=itn_irs.irs ? irs_decay(t-itn_irs.irs_time) : 0.0;

	double p_fed_irs=1.0;
	double p_fed_itn_irs=1.0;
	if(itn_irs_2016){
		const double rn=itn_rn*r_itn;
		const double rnw=itn_rnw*r_itn;
		const double dnw=itn_dnw*r_itn;
		const double ri=irs_ri*r_irs;
		const double diw=irs_diw*r_irs;
		const double riw=irs_riw*r_irs;
		p_fed_irs=(1-ri)*(1-diw-riw);
		p_fed_itn_irs=p_fed_irs*(1-rn)*(1-rnw-dnw);
	}
	else{
		const double r1=(itn_irs.itn ? itn_repel_min + (itn_repel-itn_repel_min)*r_itn : 0);
		const double d1=itn_kill*r_itn;

		double r2;
		if(irs_ellie && itn_irs.irs){
			const double ti=t-itn_irs.irs_time;
			const double det_hut_irs = 1 / (1 + exp(-(irs_decay_det1 + irs_decay_det2*ti)));
			const double succ_hut_irs = (irs_succ_k0 ? irs_k0 : 1.0) / (1 + exp(-(irs_decay_succ1 + irs_decay_succ2*ti)));
			const double mort_hut_irs = 1 / (1 + exp(-(irs_decay_mort1 + irs_decay_mort2*ti)));
			const double rep_hut_irs = 1 -  succ_hut_irs -  mort_hut_irs;

			const double kp_irs = (1 - det_hut_irs)*succ_hut_irs;
			const double jp_irs = (1 - det_hut_irs)*rep_hut_irs + det_hut_irs;
			const double lp_irs = (1 - det_hut_irs)*mort_hut_irs;
			r2=(1 - kp_irs/irs_k0)*(jp_irs/(lp_irs + jp_irs));
		}
		else if(itn_irs.irs){
			r2=irs_repel*r_irs;
		}
		else
			r2=0.0;
		p_fed_irs=(1-r2);
		p_fed_itn_irs=p_fed_irs*(1-d1-r1);
	}
	return 1-Q_in + Q_bed*p_fed_itn_irs + (Q_in-Q_bed)*p_fed_irs;
}

double Mosq_pop::seasonal_curve(const double t) const {
	return seasonal_table[modulo(static_cast<int>(floor(t)), 365)];
}

double Mosq_pop::find_eir(const double t){
	double eir = (Q0*f*I/(omega*w_bar*N));
	if(growing_pop){
		const int y=min(int(pop_r.size()-1), max(0, int(floor(t/dy + year0-init_pop_year))));
		const double q=t/dy-floor(t/dy);
		const double r=y==pop_r.size()-1 ? pop_r[y] : q*pop_r[y+1]+(1.0-q)*pop_r[y];
		eir /= r;
	}
	return eir;
}

void Mosq_pop::update(const double dt){
	const double t=village->simulation->event_manager.time_now();

	FOIM=Q*f*C/D;
	const double lagged_FOIM=lag_FOIM((FOIM+FOIM_end)/2);
	const double lagged_eir=lag_eir(find_eir(t));
	const double eir_dt = max(0.0, lagged_eir)*dt;
		
	if(eir_dt>0){
		for(int j=0; j<village->village_het.size(); j++){
			Village_het* vh=&(village->village_het[j]);
			const int nh=vh->size();
			const int num_inf=rpoisson(eir_dt*village->village_parms()->het_exp[j]*vh->size(), ran);
			for(int i=0; i<num_inf; i++){
				Human* h=(*vh)[ran_int(nh, ran)];
				h->possible_inf_bite((village->any_itn_irs ? prob_bite(h->itn_irs, t) : 1), t);	
			}
		}
	}
	FOIM_end = Q*f*C/D;
	
	//********************************************************************************
	const double mean_births=dt*(feedback ? 0.5*PL/dp : rel_M*mu0*M0*seasonal_curve(t));
	const double lambda=max(0.0, lagged_FOIM);
	const double pd=dt*mu;
	const double pe=lambda*dt;

	const double births=	(stochastic ? rpoisson(mean_births, ran) : mean_births);
	const double deaths_S=	(stochastic ? rbinom(pd, S, ran) : pd*S);
	const double inc_E=		(stochastic ? rbinom(pe/(1-pd), S-deaths_S, ran) : pe*S);
	const double lagged_inc_E=lag_inc_E(inc_E);
	const double inc_I=		(stochastic ? rbinom(PM, lagged_inc_E, ran) : PM*lagged_inc_E);
	const double deaths_I=	(stochastic ? rbinom(pd, I, ran) : pd*I);

	S += births - inc_E	- deaths_S;
	I += inc_I			- deaths_I;
	E += inc_E - inc_I - pd*E;

	S=max(S, 0.0);
	I=max(I, 0.0);
	E=max(E, 0.0);
	if (stochastic)
		E = floor(E+0.5);
	
	//********************************************************************************
	if(feedback){
		const double M=S + I + E; 

		const double expm4 = .01831564;
		const double larval_factor = (larvicide ?
			(larvi_min + (1.0-larvi_min)*(1+expm4)/(1+expm4*exp(-deprec_param*(t-larval_time))))
			: 1.0);
		const double K=K0*seasonal_curve(t)*larval_factor;
		find_beta();
		const double max_mu=5;
		const double change_EL=beta*M - min(max_mu, mue*(1+(EL+LL)/K))*EL-EL/de;
		const double change_LL=EL/de - min(max_mu, mul*(1+gamma*(EL+LL)/K))*LL - LL/dl;
		const double change_PL=LL/dl - mup*PL - PL/dp;

		EL += dt*change_EL;
		LL += dt*change_LL;
		PL += dt*change_PL;
		EL=max(EL, 0.0);
		LL=max(LL, 0.0);
		PL=max(PL, 0.0);
	}

	//********************************************************************************
}

void Mosq_pop::update_foim(const double sum_fail, const double sum_repeat, const double sum_rel_biting_rate){
	w_bar=1-Q0*sum_fail/sum_rel_biting_rate;
	z_bar=Q0*sum_repeat/sum_rel_biting_rate;		
	D=omega*village->total_pop()*(1-sum_fail/sum_rel_biting_rate);
	Q=1-(1-Q0)/w_bar;
	f=1/(tau1/(1-z_bar)+tau2);
	av=Q0*f;
	mu=-f*log(p1*w_bar*p2/(1-z_bar*p1));
	PM=exp(-mu*latmosq);
}

//***************************************************************************
//***************************************************************************
