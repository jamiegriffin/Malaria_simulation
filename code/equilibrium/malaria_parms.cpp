 //***************************************************************
//***************************************************************


void Parms::calc_dydt(const double t, Malaria_state &y, Malaria_state &dydt){
	if(run_type<2)
		set_to_zero(dydt);
	else
		for(int i=0; i<6*mosq_pop->size(); i++)
			dydt[i]=0;

	if(run_type<2)
		update_parms(y, t); 
	update_mosq(y, dydt, t);
	if(run_type<2){
		update_human(y, dydt, t);
		ageing(y, dydt);
	}
}

//***************************************************************
//***************************************************************

Parms::Parms() : itn_coverage(0), irs_coverage(0), run_type(0), 
	leave_itn_on(false), prev_age0(2), prev_age1(10) {}


void Parms::change_itn_coverage(const double p10){

	const double p1=min(max(p10, 1E-12), 1.0-1E-12);
	const double p0=itn_coverage;
	itn_coverage=p1;
	if(run_type<2){
		const int num_prev_states=state.human_state[0][0].num_prev_states();
		const int num_states=state.human_state[0][0].num_states();
		const bool increase=p1>p0;
		for(int j=0; j<num_het; j++){
			Human_state& h0=state.human_state[0][j];
			Human_state& h1=state.human_state[1][j];
			for(int k=0; k<num_prev_states; k++){
				double* s0=h0.state(k);
				double* s1=h1.state(k);
				for(int i=0; i<num_age; i++){
					const double x=increase ? s0[i]*(p1-p0) : s1[i]*(p1-p0);
					s0[i] -= x;
					s1[i] += x;
				}
			}

			for(int k=num_prev_states; k<num_states; k++){
				double* s0=h0.state(k);
				double* s1=h1.state(k);
				for(int i=0; i<num_age; i++)
					increase ? s1[i]=(s1[i]*p0+s0[i]*(p1-p0))/p1 : s0[i]=(s0[i]*(1-p0)+s1[i]*(p0-p1))/(1-p1);			
			}
		}

		for(int i=0; i<lag_mean_c[0].history.size(); i++){
			const double x=increase ? lag_mean_c[0].history[i].second*(p1-p0) : 
									lag_mean_c[1].history[i].second*(p1-p0);
			lag_mean_c[0].history[i].second -= x;
			lag_mean_c[1].history[i].second += x;
		}
		const double x=increase ? mean_c[0]*(p1-p0) : mean_c[1]*(p1-p0);
		mean_c[0] -= x;
		mean_c[1] += x;
	}

}
	
void Parms::add_itn(const double coverage){
	if(!itn_irs_det.itn){
		leave_itn_on=from_map_bool("itn_leave");
		leave_itn_time=dy*from_map("itn_leave_dur");
		state.current_size=state.total_size;
		lag_EIR[1]=lag_EIR[0];
		mean_c[1]=mean_c[0];
		lag_mean_c[1]=lag_mean_c[0];
	}

	change_itn_coverage(coverage);
	itn_irs_det.itn=true;
	itn_irs_det.itn_time=ode.time();
	update_coverage_time=itn_irs_det.itn_time;
}


void Parms::add_irs(const double coverage){
	if(run_type<2)
		error_crit("Can't introduce IRS when running the deterministic human model");
	const double t=ode.time();

	if(itn_irs_det.irs && irs_coverage>coverage+1E-3){
		double decay=0.0;
		int nmm=0;
		const int nm=mosq_pop->size();
		for(int m=0; m<nm; m++){
			if((*mosq_pop)[m].M0>0){
				nmm++;
				decay += (*mosq_pop)[m].irs_decay(t-itn_irs_det.irs_time);
			}
		}
		decay /= nmm;
		irs_coverage=coverage + (irs_coverage-coverage)*decay;
	}
	else
		irs_coverage=coverage;

	itn_irs_det.irs=true;
	itn_irs_det.irs_time=t;
}

void Parms::setup(const int na, const int nh, const vector<vector<double> > &v, const double birth_rate0){
	num_age=na;
	num_het=nh;
	prev_age0=from_map("prev_age0");
	prev_age1=from_map("prev_age1");
	input_full();
	init_demography(v, birth_rate0);
	reset();
	state.setup(*this);
}

void Parms::reset(){
	leave_itn_on=from_map_bool("itn_leave");
	init_mosq();
	init_infection();
}

void Parms::input_full(){

	ft=-1E20;
	cT=-1E20;

	inf_immunity =	from_map_int("inf_immunity", 0, 6);
	clin_immunity = from_map_int("clin_immunity", 0, 6);
	sev_immunity = from_map_int("sev_immunity", 0, 7);
	det_immunity = from_map_int("det_immunity", 0, 6);
	
	separate_sev=from_map_bool("separate_sev", true);

	eta=from_map("eta");
	rho=from_map("rho", 0, 1);
	a0=from_map("a0");
	
	dur_E=from_map("dur_E");
	dur_D=from_map("dur_D");
	dur_I =	from_map("dur_A")+dur_D;
	dur_U=from_map("dur_U");

	bh=from_map("bh", 0, 1);
	bmin=from_map("bmin", 0, 1);
	db=from_map("db");
	kb=from_map("kb");
	IB0=from_map("IB0");
	ub=from_map("ub");
	gammab=from_map("gammab");
	fb0=from_map("fb0");
	ab0=from_map("ab0");
	
	dmin=from_map("dmin", 0, 1);
	dd=from_map("dd");
	kd=from_map("kd");
	ID0=from_map("ID0");
	ud=from_map("ud");
	gammad=from_map("gammad");
	fd0=from_map("fd0");
	ad0=from_map("ad0");
	
	phi0=from_map("phi0", 0, 1);
	phi1=from_map("phi1", 0, 1);
	dc=from_map("dc");
	kc=from_map("kc");
	IC0=from_map("IC0");
	uc=from_map("uc");
	gammac=from_map("gammac");
	fc0=from_map("fc0");
	ac0=from_map("ac0");
	
	P_IC_M=from_map("P_IC_M", 0, 1);
	dm=from_map("dm");
	
	theta0=from_map("theta0", 0, 1);
	theta1=from_map("theta1", 0, 1);
	dv=from_map("dv");
	kv=from_map("kv");
	IV0=from_map("IV0");
	uv=from_map("uv");
	gammav=from_map("gammav");
	fv0=from_map("fv0");
	av0=from_map("av0");
	P_IV_M=from_map("P_IV_M", 0, 1);
	dvm=from_map("dvm");
	tau_v=in_map("tau_v") ? from_map("tau_v") : 1.0;
	
	sigma2=from_map("sigma2");
	
	latgam=from_map("latgam");

	cD=from_map("cD", 0, 1);
	cU=from_map("cU", 0, 1);
	beta_pcr=from_map("beta_pcr");
	
	gamma_inf=from_map("gamma_inf");
	alpha_pcr=from_map("alpha_pcr");
}

void Parms::update_history(Malaria_state &y, const double current_time){
	const int nh=(itn_coverage>0 ? 2 : 1);
	const int num_prev_states=state.human_state[0][0].num_prev_states();
	if(run_type<2){
		for(int j=0; j<num_het; j++){
			for(int i=0; i<num_age; i++){
				double Z=0;
				for(int k=0; k<nh; k++){
					Human_state& h=y.human_state[k][j];
					for(int l=0; l<num_prev_states; l++)
						Z += h.state(l)[i];
				}
				const double r=N_pop[i]*het_wt[j]/Z;
				for(int k=0; k<nh; k++){
					Human_state& h=y.human_state[k][j];
					for(int l=0; l<num_prev_states; l++)
						h.state(l)[i] *= r;
				}
			}
		}
	}

	update_foim(true, current_time);

	if(run_type<2){
		lag_prev.update(current_time, par_prev(y, prev_age0*dy, prev_age1*dy));
		find_eir(y, current_time);
		for(int k=0; k<nh; k++)
			lag_EIR[k].update(current_time, EIR[k]);
		update_parms(y, current_time);	
		update_sum_H(y);
	}
	const int nm=mosq_pop->size();
	for(int m=0; m<nm; m++)
		if((*mosq_pop)[m].M0>0)
			(*mosq_pop)[m].update_history(current_time);

	if(itn_coverage>0 && leave_itn_on && update_coverage_time<current_time-1){ 
		change_itn_coverage(exp(-(current_time-update_coverage_time)/leave_itn_time)*itn_coverage);
		update_coverage_time=current_time;	
	}
}
void Parms::clear_history(){
	for(int k=0; k<lag_mean_c.size(); k++){
		lag_EIR[k].clear();
		lag_mean_c[k].clear();
	}
	lag_prev.clear();

	for(int m=0; m<mosq_pop->size(); m++)
		if((*mosq_pop)[m].M0>0)
			(*mosq_pop)[m].lag_infv.clear();
}
void Parms::rotate_history(const double t){
	for(int k=0; k<lag_mean_c.size(); k++){
		lag_EIR[k].rotate(t);
		lag_mean_c[k].rotate(t);
	}
	lag_prev.rotate(t);

	for(int m=0; m<mosq_pop->size(); m++)
		if((*mosq_pop)[m].M0>0)
			(*mosq_pop)[m].lag_infv.rotate(t);
}
void Parms::copy_history(const Parms &parms){
	for(int k=0; k<lag_mean_c.size(); k++){
		lag_EIR[k]=parms.lag_EIR[k];
		lag_mean_c[k]=parms.lag_mean_c[k];
	}
	lag_prev=parms.lag_prev;

	for(int m=0; m<mosq_pop->size(); m++)
		if((*mosq_pop)[m].M0>0)
			(*mosq_pop)[m].lag_infv=(*parms.mosq_pop)[m].lag_infv;
}

void Parms::find_eir(const Malaria_state &y, const double t){
	for(int k=0; k<(itn_coverage>0 ? 2 : 1); k++){
		EIR[k]=0;
		for(int m=0; m<mosq_pop->size(); m++)
			if((*mosq_pop)[m].M0>0)
				EIR[k] += (*mosq_pop)[m].find_eir(t)*(k==0 ? 1 : (*mosq_pop)[m].prob_bite(itn_irs_det, t));
	}
}

void Parms::update_parms(const Malaria_state &y, const double t){
	for(int k=0; k<(itn_coverage>0 ? 2 : 1); k++){
		for(int j=0; j<num_het; j++){
			const Human_state &h=y.human_state[k][j];
			IC_20[k][j]=h.IC_A[age_20_l] + age_20_factor*(h.IC_A[age_20_u]-h.IC_A[age_20_l]);
			IV_20[k][j]=h.IV[age_20_l] + age_20_factor*(h.IV[age_20_u]-h.IV[age_20_l]);
		}
	}
}
void Parms::update_mosq(const Malaria_state &y,  Malaria_state &dydt, const double t){
	const int nm=mosq_pop->size();
	for(int m=0; m<nm; m++)
		if((*mosq_pop)[m].M0>0)
			(*mosq_pop)[m].update_det(y.mosq_state[m], dydt.mosq_state[m], t);
}
void Parms::update_foim(const bool update_lag, const double t){
	const int nh=(itn_coverage>0 ? 2 : 1);
	const int nm=mosq_pop->size();

	if(run_type<2){
		find_mean_c(state);
		if(update_lag)
			for(int k=0; k<nh; k++)
				lag_mean_c[k].update(t, mean_c[k]);
	}

	for(int m=0; m<nm; m++){
		if((*mosq_pop)[m].M0>0){
			const bool itn=itn_irs_det.itn;
			const bool irs=itn_irs_det.irs;
			double sum_fail, sum_repeat;
			if(!itn && !irs){
				if(run_type<2)
					(*mosq_pop)[m].C=(*mosq_pop)[m].N*lag_mean_c[0](t);
				sum_fail=0.0;
				sum_repeat=0.0;
			}
			else if(!irs){
				const double prob_fail=(*mosq_pop)[m].prob_fail(itn_irs_det, t);
				if(run_type<2)
					(*mosq_pop)[m].C=(*mosq_pop)[m].N*(lag_mean_c[0](t) + lag_mean_c[1](t)*(1-prob_fail));
				sum_fail=itn_coverage*prob_fail;
				sum_repeat=itn_coverage*(*mosq_pop)[m].prob_repeat(itn_irs_det, t);
			}
			else if(!itn){
				sum_fail=irs_coverage*(*mosq_pop)[m].prob_fail(itn_irs_det, t);
				sum_repeat=irs_coverage*(*mosq_pop)[m].prob_repeat(itn_irs_det, t);
			}
			else{
				itn_irs_det.irs=false;
				const double pf10=(*mosq_pop)[m].prob_fail(itn_irs_det, t);
				const double pr10=(*mosq_pop)[m].prob_repeat(itn_irs_det, t);
				itn_irs_det.irs=true;
				itn_irs_det.itn=false;
				const double pf01=(*mosq_pop)[m].prob_fail(itn_irs_det, t);
				const double pr01=(*mosq_pop)[m].prob_repeat(itn_irs_det, t);
				itn_irs_det.itn=true;
				const double pf11=(*mosq_pop)[m].prob_fail(itn_irs_det, t);
				const double pr11=(*mosq_pop)[m].prob_repeat(itn_irs_det, t);

				const double p0=itn_coverage;
				const double p1=irs_coverage;
				
				sum_fail=p0*(1-p1)*pf10 + (1-p0)*p1*pf01 + p0*p1*pf11;
				sum_repeat=p0*(1-p1)*pr10 + (1-p0)*p1*pr01 + p0*p1*pr11;
			}
			(*mosq_pop)[m].update_foim(sum_fail, sum_repeat, 1.0);
		}
	}
}
void Parms::update_human(const Malaria_state &y, Malaria_state &dydt, const double t){
	set_to_zero(clin_inc);
	set_to_zero(inf_inc);
	set_to_zero(sev_inc);
	for(int k=0; k<(itn_coverage>0 ? 2 : 1); k++){
		const double lagged_eir=lag_EIR[k](t);
		for(int j=0; j<num_het; j++){
			const double eird_max=lagged_eir*rel_foi[j];
			update_human(y.human_state[k][j], dydt.human_state[k][j], j, k, t, eird_max);
		}
	}
}

void Parms::update_human(const Human_state &y, Human_state &dydt, const int j, const int k, const double t, const double eird_max){
	for(int i=0; i<N_pop.size(); i++){
		const double FOIhi=find_foi(eird_max, i, y.IB[i]);
		const double phi=find_phi(i, j, k, y.IC_A[i]);
		const double AT = y.A1[i]+y.U[i];
		const double Y=y.S[i]+AT;
		const double Fp=phi*FOIhi;
		const double inc_symp=Fp*Y;
		inf_inc[i] += FOIhi*Y/N_pop[i];
		clin_inc[i] += inc_symp/N_pop[i];
		sev_inc[i] += FOIhi*find_theta(i, j, k, y.IV[i])*Y/N_pop[i]*(separate_sev ? 1.0 : phi);

	//****************************************************
	//****************************************************

		const double infection		=	FOIhi*y.S[i];
		const double recoveryT		=	rec_T*y.T[i];
		const double end_P			=	rec_P*y.P[i];
		const double recoveryD		=	rec_D*y.D[i];
		const double recoveryA1		=	rec_A0*y.A1[i];
		const double recoveryU		=	rec_U*y.U[i];
		dydt.S[i]	+=	-infection + end_P + recoveryU;
		dydt.T[i]	+=	ft*inc_symp -  recoveryT;
		dydt.D[i]	+=	(1-ft)*inc_symp -  recoveryD;
		dydt.A1[i]	+=	FOIhi*(1-phi)*(y.U[i]+y.S[i]) - phi*FOIhi*y.A1[i] + recoveryD-recoveryA1;
		dydt.U[i]	+=	-FOIhi*y.U[i] + recoveryA1 - recoveryU;	
		dydt.P[i]	+=	recoveryT - end_P;
		
		if(inf_immunity>0)
			dydt.IB[i] += increase_IB(eird_max*FOIh_age[i])	- dec_db*y.IB[i];
		if(clin_immunity>0)
			dydt.IC_A[i] += increase_IC(FOIhi)	- dec_dc*y.IC_A[i];
		if(sev_immunity>0)
			dydt.IV[i] += increase_IV(FOIhi)	- dec_dv*y.IV[i];
		if(det_immunity>0)
			dydt.ID[i] += increase_ID(FOIhi)	- dec_dd*y.ID[i];
	}
}

void Parms::init_mosq(){
	lag_EIR.resize(2);
	lag_mean_c.resize(2);
	for(int k=0; k<lag_mean_c.size(); k++){
		lag_EIR[k].setup(dur_E);
		lag_mean_c[k].setup(latgam);
	}
	EIR.assign(2, 0);
	mean_c.assign(2, 0);
}
void Parms::init_infection(){

	lag_prev.setup(1.1*dy);
	const int n_table=10000;
	b_table.setup(Func_hill(bh, bmin, IB0, kb),			n_table, 0.00001*IB0, 100*IB0);
	c_table.setup(Func_hill(phi0, phi1, IC0, kc),		n_table, 0.00001*IC0, 100*IC0);
	d_table.setup(Func_hill(1.0, dmin, ID0, kd),		n_table, 0.00001*ID0, 100*ID0);
	v_table.setup(Func_hill(theta0, theta1, IV0, kv),	n_table, 0.00001*IV0, 100*IV0);

	if(rel_foi.size() != num_het){
		het_x.resize(num_het);
		het_wt.resize(num_het);

		if(num_het>9){
			het_wt.assign(het_wt.size(), 1.0/het_wt.size());
			for(int j=0; j<het_x.size(); j++)
				het_x[j]=invnorm((j+0.5)/het_x.size());
		}
		else
			gauher_norm(het_x, het_wt);
		rel_foi.resize(num_het);
	}
	if(rel_foi.size()==1)
		rel_foi[0]=1;
	else
		for(int j=0; j<rel_foi.size(); j++)
			rel_foi[j]=exp(-sigma2/2+sqrt(sigma2)*het_x[j]);
	rec_A0=1/(dur_I-dur_D);
	rec_D=1/dur_D;
	rec_U=1/dur_U;
	rec_T=1/dur_T;
	rec_P=1/(dur_P-dur_T);

	IC_20.assign(2, vector<double>(num_het, 0));
	IV_20.assign(2, vector<double>(num_het, 0));
	
	if(IC_M_age.size() != N_pop.size())
		IC_M_age.resize(N_pop.size());
	if(IV_M_age.size() != N_pop.size())
		IV_M_age.resize(N_pop.size());

	fB.assign(N_pop.size(), 1);
	fC.assign(N_pop.size(), 1);
	fD.assign(N_pop.size(), 1);
	fV.assign(N_pop.size(), 1);
	if(inf_immunity>0){
		dec_db=1/db;
		if(inf_immunity==5 || inf_immunity==6)
			for(int i=0; i<N_pop.size(); i++)
				fB[i]=1-(1-fb0)/(1+pow((i==N_pop.size()-1 ? age[i] : 0.5*(age[i]+age[i+1]))/ab0, gammab));
	}
	dec_dc=1/dc;
	for(int i=0; i<N_pop.size(); i++)
		IC_M_age[i]=i==N_pop.size()-1 ? 0.0 : P_IC_M*dm/(age[i+1]-age[i])*(exp(-age[i]/dm)-exp(-age[i+1]/dm));

	if(clin_immunity==5 || clin_immunity==6)
		for(int i=0; i<N_pop.size(); i++)
			fC[i]=1-(1-fc0)/(1+pow((i==N_pop.size()-1 ? age[i] : 0.5*(age[i]+age[i+1]))/ac0, gammac));



	if(det_immunity>0){
		dec_dd=1/dd;
		if(det_immunity==5 || det_immunity==6)
			for(int i=0; i<N_pop.size(); i++)
				fD[i]=1-(1-fd0)/(1+pow((i==N_pop.size()-1 ? age[i] : 0.5*(age[i]+age[i+1]))/ad0, gammad));
	}
	if(sev_immunity>0){
		dec_dv=1/dv;
	for(int i=0; i<N_pop.size(); i++)
		IV_M_age[i]=i==N_pop.size()-1 ? 0.0 : P_IV_M*dvm/(age[i+1]-age[i])*(exp(-age[i]/dvm)-exp(-age[i+1]/dvm));
		if(sev_immunity==5 || sev_immunity==6 || sev_immunity==7)
			for(int i=0; i<N_pop.size(); i++)
				fV[i]=1-(1-fv0)/(1+pow((i==N_pop.size()-1 ? age[i] : 0.5*(age[i]+age[i+1]))/av0, gammav));
	}
	
	if(FOIh_age.size() != N_pop.size()){
		FOIh_age.resize(N_pop.size());
		clin_inc.resize(N_pop.size());
		inf_inc.resize(N_pop.size());
		sev_inc.resize(N_pop.size());
	}
	for(int i=0; i<N_pop.size(); i++)
		FOIh_age[i]=1-rho*exp(-(i==N_pop.size()-1 ? age[i] : 0.5*(age[i]+age[i+1]))/a0);

	vector<double> N_foi(N_pop.size());
	for(int i=0; i<N_pop.size(); i++)
		N_foi[i]=N_pop[i]*FOIh_age[i];
	omega=sum(N_foi);

	clin_inc.assign(N_pop.size(), 0);
	inf_inc.assign(N_pop.size(), 0);
	sev_inc.assign(N_pop.size(), 0);

	fB1=fB;
	fC1=fC;
	fD1=fD;
	fV1=fV;
	for(int i=0; i<N_pop.size(); i++){
		fB1[i]=pow(fB[i], 1.0/kb);
		fC1[i]=pow(fC[i], 1.0/kc);
		fD1[i]=pow(fD[i], 1.0/kd);
		fV1[i]=pow(fV[i], 1.0/kv);
	}

}

double Parms::increase_IB(const double eir1){
	return inf_immunity==0 ? 0 :	(inf_immunity==1  ? 1 : ((inf_immunity==4 || inf_immunity==6) ? eir1/(eir1*ub+1) : eir1));
}
double Parms::increase_IC(const double FOIh1){
	return clin_immunity==0 ? 0 :	(clin_immunity==1  ? 1 : ((clin_immunity==4 || clin_immunity==6) ? FOIh1/(FOIh1*uc+1) : FOIh1));
}
double Parms::increase_IV(const double FOIh1){
	return sev_immunity==0 ? 0 :	(sev_immunity==1  ? 1 : ((sev_immunity==4 || sev_immunity==6) ? FOIh1/(FOIh1*uv+1) : ((sev_immunity==3 || sev_immunity==7) ?  pow(FOIh1, tau_v) : FOIh1)));
}
double Parms::increase_ID(const double FOIh1){
	return det_immunity==0 ? 0 :	(det_immunity==1  ? 1 : ((det_immunity==4 || det_immunity==6) ? FOIh1/(FOIh1*ud+1) : FOIh1));
}

double Parms::find_phi(const int i, const int j, const int k, const double ICA){
	return clin_immunity==0 ? phi0 : c_table((IC_M_age[i]*IC_20[k][j]+ICA)*fC1[i]);
}
double Parms::find_theta(const int i, const int j, const int k, const double IV1){
	return sev_immunity==0 ? theta0 : v_table((IV_M_age[i]*IV_20[k][j]+IV1)*fV1[i]);
}
double Parms::find_b(const int i, const double IB){
	return (inf_immunity==0 || IB==0) ? bh : b_table(IB*fB1[i]);
}
double Parms::prob_det(const int i, const double ID) const{
	return (det_immunity==0 || ID==0) ? 1.0 : d_table(ID*fD1[i]);
}
double Parms::find_foi(const double eird_max, const int i, const double IB){
	return max(0.0, eird_max*FOIh_age[i]*find_b(i, IB));
}

void Parms::find_mean_c(const Malaria_state &y){
	const int N=N_pop.size();
	for(int k=0; k<(itn_coverage>0 ? 2 : 1); k++){
		mean_c[k]=0;
		for(int j=0; j<num_het; j++){
			const Human_state& h=y.human_state[k][j];
			double cj=0.0;
			for(int i=0; i<N; i++){
				const double p_det=prob_det(i, h.ID[i]);
				const double cji=cT*h.T[i]+cD*h.D[i]+(cU + (cD-cU)*pow(p_det, gamma_inf))*h.A1[i]+cU*h.U[i];
				cj += cji*FOIh_age[i];
			}
			mean_c[k] += cj*rel_foi[j];
		}
	}
}

//***************************************************************

void Parms::human_equilibrium(Human_state &y, const double eird_max, const int j){

	if(eird_max<1E-11){
		set_to_zero(y);
		for(int i=0; i<N_pop.size(); i++)
			y.S[i]=N_pop[i]*het_wt[j];
	}
	else{
 		for(int i=0; i<N_pop.size(); i++){
			const double eiri=eird_max*FOIh_age[i];
			if(inf_immunity>0)
				y.IB[i]=((i==0 ? 0 : y.IB[i-1])+ increase_IB(eiri)*x_I[i])/(1+x_I[i]/db);

			const double FOIhi=find_foi(eird_max, i, y.IB[i]);
			if(clin_immunity>0)
				y.IC_A[i]=	((i==0 ? 0 : y.IC_A[i-1])+ increase_IC(FOIhi)*x_I[i])/(1+x_I[i]/dc);
			if(sev_immunity>0)
				y.IV[i]=	((i==0 ? 0 : y.IV[i-1])+ increase_IV(FOIhi)*x_I[i])/(1+x_I[i]/dv);
			if(det_immunity>0)
				y.ID[i]=	((i==0 ? 0 : y.ID[i-1])+ increase_ID(FOIhi)*x_I[i])/(1+x_I[i]/dd);
		}
		if(clin_immunity>0)
			IC_20[0][j]=y.IC_A[age_20_l]+age_20_factor*(y.IC_A[age_20_u]-y.IC_A[age_20_l]);
		if(sev_immunity>0)
			IV_20[0][j]=y.IV[age_20_l]+age_20_factor*(y.IV[age_20_u]-y.IV[age_20_l]);

		for(int i=0; i<N_pop.size(); i++){
			const double FOIhi=find_foi(eird_max, i, y.IB[i]);
			const double phi=find_phi(i, j, 0, y.IC_A[i]);

			const double delta= (i==0 ? eta : age_rate[i-1]);
			const double gamma= eta+ (i==N_pop.size()-1 ? 0 : age_rate[i]);
			
			const double betaS=FOIhi + gamma;
			const double betaT=rec_T + gamma;
			const double betaP=rec_P + gamma;
	
			const double FOIi_T=ft*phi*FOIhi;
			const double FOIi_D=(1-ft)*phi*FOIhi;
			const double FOIi_asymp=(1-phi)*FOIhi;
						
			const double S_prev=	(i==0 ? het_wt[j] : y.S[i-1]);
			const double P_prev=	(i==0 ? 0 : y.P[i-1]);
			const double T_prev=	(i==0 ? 0 : y.T[i-1]);
			const double D_prev=	(i==0 ? 0 : y.D[i-1]);
			const double A1_prev=	(i==0 ? 0 : y.A1[i-1]);
			
			const double aT=FOIi_T/betaT;
			const double bT=delta*T_prev/betaT;
			const double aP=rec_T*aT/betaP;
			const double bP=(rec_T*bT+delta*P_prev)/betaP;
			
			double Y;
			const double U_prev=	(i==0 ? 0 : y.U[i-1]);
			const double Y_prev=S_prev+A1_prev+U_prev;

			const double betaD=rec_D + gamma;
			const double betaA1=FOIhi*phi + rec_A0 + gamma;
			const double betaU= FOIhi + rec_U + gamma;
							
			const double aD=FOIi_D/betaD;
			const double bD=delta*D_prev/betaD;
			Y=(N_pop[i]*het_wt[j]-(bT+bD+bP))/(1+aT+aD+aP);
			y.T[i]=(aT*Y+bT);
			y.D[i]=(aD*Y+bD);
			y.P[i]=(aP*Y+bP);

			y.A1[i]=(delta*A1_prev+FOIi_asymp*Y+rec_D*y.D[i])/(betaA1+FOIi_asymp);
			y.U[i]=(rec_A0*y.A1[i]+delta*U_prev)/betaU;
			y.S[i]=Y - y.A1[i] -y.U[i];
			const double inf_inc_ij=FOIhi*Y/N_pop[i];
			clin_inc[i] += inf_inc_ij*phi; 	
			inf_inc[i] += inf_inc_ij;
			sev_inc[i] += inf_inc_ij*find_theta(i, j, 0, y.IV[i])*(separate_sev ? 1.0 : phi);	
		}
	}
}

void Parms::mosq_equilibrium(Malaria_state &y){
	find_mean_c(y);
	for(int m=0; m<mosq_pop->size(); m++)
		(*mosq_pop)[m].find_equilibrium_det(y.mosq_state[m], mean_c[0]/omega);
}

void Parms::equilibrium_foi(Malaria_state &y, const double eird_max){
	set_to_zero(y);
	clin_inc.assign(clin_inc.size(), 0);
	inf_inc.assign(inf_inc.size(), 0);
	sev_inc.assign(sev_inc.size(), 0);
	
	for(int j=0; j<rel_foi.size(); j++)
		human_equilibrium(y.human_state[0][j], rel_foi[j]*eird_max, j);

	mosq_equilibrium(y);
}
void Parms::update_sum_H(const Malaria_state &y){
	Sh=0;
	Th=0;
	Dh=0;
	A1h=0;
	Ph=0;
	Uh=0;
	for(int k=0; k<(itn_coverage>0 ? 2 : 1); k++){
		for(int j=0; j<num_het; j++){
			const Human_state &h=y.human_state[k][j];
			Sh += sum(h.S, num_age);
			Th += sum(h.T, num_age);
			Dh += sum(h.D, num_age);
			A1h += sum(h.A1, num_age);
			Ph += sum(h.P, num_age);
			Uh += sum(h.U, num_age);
			IC_20[k][j]=h.IC_A[age_20_l]+age_20_factor*(h.IC_A[age_20_u]-h.IC_A[age_20_l]);
			IV_20[k][j]=h.IV[age_20_l]+age_20_factor*(h.IV[age_20_u]-h.IV[age_20_l]);
		}
	}
	H_sum=Sh+Th+Dh+A1h+Uh+Ph;
}

void Parms::update_from_equilibrium(Malaria_state &y){
	update_sum_H(y);
	update_history(y, ode.time());
}
//******************************************************************************
//******************************************************************************

void Parms::find_equilibrium(Malaria_state &y){
	set_to_zero(y);
	double mv0=0;
	for(int m=0; m<mosq_pop->size(); m++)
		mv0 += (*mosq_pop)[m].M0*(*mosq_pop)[m].rel_M/(*mosq_pop)[m].N;
	const double tol=1E-6;
	double diff=1;
	double eird_max=5*mv0/365.0;
	double eird_max_old;
	while(diff > tol){	
		equilibrium_foi(y, eird_max);
		eird_max_old=eird_max;
		find_eir(y);
		eird_max=EIR[0];
		diff=abs(eird_max-eird_max_old)/eird_max;
	}
	update_from_equilibrium(y);
	EIR[0]=eird_max;
}

void Parms::find_equilibrium_known_eir(Malaria_state &y, const double eiry0){

	set_to_zero(y);
	const double eird_max=eiry0/dy;
	clin_inc.assign(clin_inc.size(), 0);
	inf_inc.assign(inf_inc.size(), 0);
	sev_inc.assign(sev_inc.size(), 0);
	for(int j=0; j<rel_foi.size(); j++)
		human_equilibrium(y.human_state[0][j], rel_foi[j]*eird_max, j);

	find_mean_c(y);	

	update_sum_H(y);
	update_history(y, 0);
	EIR[0]=eird_max;
}

//***************************************************************
//***************************************************************
double Parms::par(const Human_state &y, const int i, const bool pcr) const{
	const double p_det=prob_det(i, y.ID[i]);
	return y.T[i] + y.D[i] + (pcr ? pow(p_det, alpha_pcr)*y.A1[i]  + pow(p_det, beta_pcr)*y.U[i] : p_det*y.A1[i]);
}
double Parms::par_prev(const Malaria_state &y, const int i, const bool pcr) const{
	double prev=0;
	for(int k=0; k<(itn_coverage>0 ? 2 : 1); k++){
		for(int j=0; j<num_het; j++)
			prev += par(y.human_state[k][j], i, pcr)/N_pop[i];
	}
	return prev;
}
void Parms::par_prev(const Malaria_state &y, vector<double> &pred, const bool pcr) const{
	const int N=N_pop.size();
	for(int i=0; i<N; i++)
		pred[i]=par_prev(y, i, pcr);
}
double Parms::par_prev(const Malaria_state &y, const double age1, const double age2, const bool pcr) const{
	set_to_zero(temp_pred);
	const int i1=find_age_group(age1);
	const int i2=min(int(N_pop.size()), find_age_group(age2)+1);
	for(int i=i1; i<i2; i++)
		temp_pred[i]=par_prev(y, i, pcr);
	return mean_in_age_range(temp_pred, age1, age2);
}
double Parms::par_prev(const Malaria_state &y, const bool pcr) const{
	return par_prev(y, 0, 200*dy, pcr);
}
//***************************************************************
//***************************************************************