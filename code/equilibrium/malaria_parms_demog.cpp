
//***************************************************************
//***************************************************************
void Parms::ageing(const Malaria_state &y,  Malaria_state &dydt){
	const double births=birth_rate;

	const int nh=(itn_coverage>0 ? 2 : 1);
	const int ns=y.human_state[0][0].num_prev_states();
	for(int k=0; k<nh; k++){
		for(int j=0; j<num_het; j++){
			const Human_state &h=y.human_state[k][j];
			Human_state &dh=dydt.human_state[k][j];
			dh.S[0] += births*het_wt[j]*(k==1 ? itn_coverage : (1-itn_coverage));

			for(int l=0; l<ns; l++)
				ageing_and_death(h.state(l), dh.state(l));
			
			if(clin_immunity>0)
				ageing_mean(h.IC_A, dh.IC_A);
			if(inf_immunity>0)
				ageing_mean(h.IB, dh.IB);
			if(sev_immunity>0)
				ageing_mean(h.IV, dh.IV);
			if(det_immunity>0)
				ageing_mean(h.ID, dh.ID);
		}
	}
}

int Parms::find_age_group(const double age1) const{
	const double age2=age1+1E-8;
	return (age2>=age.back() ? age.size()-1 : age_guide(age2/age.back()));
}
double Parms::prop_in_age_range(const double age1, const double age2) const{
	const int N=N_pop.size();
	if(age2<age1)
		error_crit("Second age needs to be at least as big as first age in Parms::prop_in_age_range");

	const int r=find_age_group(age1);
	const int s=find_age_group(age2);
	if(r==s){
		if(r<N-1)
			return N_pop[r]*(age2-age1)/b[r];
		else
			return N_pop[r]*(exp(-age1*alpha.back())-exp(-age2*alpha.back()));
	}

	double S=0;
	for(int i=r; i<=s; i++)
		S += (i==r ? (age[r+1]-age1)/b[r]*N_pop[r] : (i<s ? N_pop[i] : 
					((s==N-1) ? N_pop[s]*(1.0-exp(-eta*age2))  : (age2-age[s])/b[s]*N_pop[s])));
	return S;
}
double Parms::mean_in_age_range(const vector<double> &mu, const double age1, const double age2) const{
	const int N=N_pop.size();
	if(mu.size() != N)
		error_crit("Input vector needs to be the same size as number of age groups in Parms::mean_in_age_range");
	if(age2<age1)
		error_crit("Second age needs to be at least as big as first age in Parms::mean_in_age_range");

	const int r=find_age_group(age1);
	const int s=find_age_group(age2);
	if(r==s)
		return mu[r];

	double Smu=0;
	double Sw=0;
	for(int i=r; i<=s; i++){
		const double w=(i==r ? (age[r+1]-age1)/b[r]*N_pop[r] : (i<s ? N_pop[i] : ((s==N-1) ? exp(-eta*age[N-1])-exp(-eta*age2)  : (age2-age[s])/b[s]*N_pop[s])));
		Smu += w*mu[i];
		Sw += w;
	}
	return Smu/Sw;
}

void Parms::ageing_and_death(const double *y, double *dydt) const{
	for(int i=0; i<num_age; i++)
		dydt[i] += (i==0 ? 0.0 : y[i-1]*age_rate[i-1]) - y[i]*(age_rate[i] + alpha[i]);
}

void Parms::ageing_mean(const double *y, double *dydt, const double y_at_birth) const{
	dydt[0] += (y_at_birth-y[0])*birth_rate/N_pop[0];
	for(int i=1; i<num_age; i++)
		dydt[i] += N_pop[i-1]*age_rate[i-1]/N_pop[i]*(y[i-1]-y[i]);
}

void Parms::init_demography(const vector<vector<double> > &v, const double birth_rate0){
	
	birth_rate=birth_rate0/365.0;

	b.resize(num_age-1);
	age.resize(num_age);
	N_pop.resize(num_age);
	int ia=0;
	if(num_age==36){
		for(int j=0; j<8; j++)
			b[ia++]=0.25*dy;
		for(int j=0; j<8; j++)
			b[ia++]=0.5*dy;
		for(int j=0; j<9; j++)
			b[ia++]=dy;
		for(int j=0; j<2; j++)
			b[ia++]=2.5*dy;
		for(int j=0; j<8; j++)
			b[ia++]=5*dy;
	}
	else if(num_age==44){
		for(int j=0; j<8; j++)
			b[ia++]=0.25*dy;
		for(int j=0; j<8; j++)
			b[ia++]=0.5*dy;
		for(int j=0; j<9; j++)
			b[ia++]=dy;
		for(int j=0; j<2; j++)
			b[ia++]=2.5*dy;
		for(int j=0; j<16; j++)
			b[ia++]=5*dy;
	}
	else if(num_age==34){
		for(int j=0; j<8; j++)
			b[ia++]=0.25*dy;
		for(int j=0; j<8; j++)
			b[ia++]=0.5*dy;
		for(int j=0; j<4; j++)
			b[ia++]=dy;
		for(int j=0; j<5; j++)
			b[ia++]=2*dy;
		for(int j=0; j<8; j++)
			b[ia++]=5*dy;
	}
	else if(num_age==41){
		for(int j=0; j<10; j++)
			b[ia++]=dy/10.0;
		for(int j=0; j<5; j++)
			b[ia++]=dy/5.0;
		for(int j=0; j<8; j++)
			b[ia++]=0.5*dy;
		for(int j=0; j<4; j++)
			b[ia++]=dy;
		for(int j=0; j<5; j++)
			b[ia++]=2*dy;
		for(int j=0; j<8; j++)
			b[ia++]=5*dy;
	}
	else if(num_age==86){
		for(int j=0; j<40; j++)
			b[ia++]=0.25*dy;
		for(int j=0; j<20; j++)
			b[ia++]=0.5*dy;
		for(int j=0; j<10; j++)
			b[ia++]=dy;
		for(int j=0; j<15; j++)
			b[ia++]=2*dy;
	}
	else if(num_age==137){
		for(int j=0; j<40; j++)
			b[ia++]=0.05*dy;
		for(int j=0; j<12; j++)
			b[ia++]=1.0/12.0*dy;
		for(int j=0; j<10; j++)
			b[ia++]=0.1*dy;
		for(int j=0; j<8; j++)
			b[ia++]=0.125*dy;
		for(int j=0; j<5; j++)
			b[ia++]=0.2*dy;
		for(int j=0; j<16; j++)
			b[ia++]=0.25*dy;
		for(int j=0; j<20; j++)
			b[ia++]=0.5*dy;
		for(int j=0; j<10; j++)
			b[ia++]=dy;
		for(int j=0; j<15; j++)
			b[ia++]=2*dy;
	}
	else if(num_age==145){
		for(int j=0; j<25; j++)
			b[ia++]=0.04*dy;
		for(int j=0; j<20; j++)
			b[ia++]=0.05*dy;
		for(int j=0; j<15; j++)
			b[ia++]=1.0/15.0*dy;
		for(int j=0; j<10; j++)
			b[ia++]=0.1*dy;
		for(int j=0; j<8; j++)
			b[ia++]=0.125*dy;
		for(int j=0; j<5; j++)
			b[ia++]=0.2*dy;
		for(int j=0; j<16; j++)
			b[ia++]=0.25*dy;
		for(int j=0; j<20; j++)
			b[ia++]=0.5*dy;
		for(int j=0; j<10; j++)
			b[ia++]=dy;
		for(int j=0; j<15; j++)
			b[ia++]=2*dy;
	}
	else if(num_age==149){
		for(int j=0; j<24; j++)
			b[ia++]=1.0/24.0*dy;
		for(int j=0; j<20; j++)
			b[ia++]=0.05*dy;
		for(int j=0; j<15; j++)
			b[ia++]=1.0/15.0*dy;
		for(int j=0; j<10; j++)
			b[ia++]=0.1*dy;
		for(int j=0; j<8; j++)
			b[ia++]=0.125*dy;
		for(int j=0; j<5; j++)
			b[ia++]=0.2*dy;
		for(int j=0; j<16; j++)
			b[ia++]=0.25*dy;
		for(int j=0; j<20; j++)
			b[ia++]=0.5*dy;
		for(int j=0; j<10; j++)
			b[ia++]=dy;
		for(int j=0; j<12; j++)
			b[ia++]=2.5*dy;
		for(int j=0; j<8; j++)
			b[ia++]=5*dy;
	}
	else if(num_age==18){
		for(int j=0; j<4; j++)
			b[ia++]=0.5*dy;
		for(int j=0; j<4; j++)
			b[ia++]=dy;
		for(int j=0; j<3; j++)
			b[ia++]=2*dy;
		for(int j=0; j<2; j++)
			b[ia++]=4*dy;
		for(int j=0; j<4; j++)
			b[ia++]=10*dy;
	}
	else if(num_age==11){
		for(int j=0; j<2; j++)
			b[ia++]=dy;
		for(int j=0; j<4; j++)
			b[ia++]=2*dy;
		for(int j=0; j<2; j++)
			b[ia++]=5*dy;
		for(int j=0; j<2; j++)
			b[ia++]=20*dy;
	}
	else
		error_crit(to_string(num_age) + " age groups not implemented");

	if(ia != num_age-1)
		error_crit("Age widths b[] implemented incorrectly with num_age = " + to_string(num_age));

	age_rate.resize(num_age);
	for(int i=0; i<b.size(); i++)
		age_rate[i]=1/b[i];
	age_rate[num_age-1]=0;

	age[0]=0;
	for(int i=1; i<num_age; i++)
		age[i]=b[i-1]+age[i-1];
	age_guide.setup(b);

	age_20_u=0;
	while(age[age_20_u]+0.5*b[age_20_u]<20*dy){
		age_20_u++;
		if(age_20_u == num_age-1)
			error_crit("vector age needs to go past age 20 in Parms");
	}
	age_20_l=age_20_u-1;
	age_20_factor=(20*dy-age[age_20_l]-0.5*b[age_20_l])*2.0/(b[age_20_l]+b[age_20_u]);

	alpha.assign(num_age, -1);
	for(int j=0; j<v[0].size(); j++){
		const int i=find_age_group(v[0][j]*365.0);
		alpha[i]=v[1][j]/365.0;
	}
	alpha[0]=v[1][0]/365.0;
	for(int i=1; i<num_age; i++){
		if(alpha[i]<0)
			alpha[i]=alpha[i-1];
	}

	N_pop[0]=birth_rate/(alpha[0]+age_rate[0]);
	for(int i=1; i<num_age; i++)
		N_pop[i]=age_rate[i-1]*N_pop[i-1]/(age_rate[i]+alpha[i]);

	for(int j=0; j<20; j++){
		const double Sn=sum(N_pop); 
		for(int i=0; i<num_age; i++)
			N_pop[i] /= Sn;
		double S=0;
		for(int i=0; i<num_age; i++)
			S += N_pop[i]*alpha[i];
		for(int i=0; i<num_age; i++)
			alpha[i] *= birth_rate/S;
		N_pop[0]=birth_rate/(alpha[0]+age_rate[0]);
		for(int i=1; i<num_age; i++)
			N_pop[i]=age_rate[i-1]*N_pop[i-1]/(age_rate[i]+alpha[i]);
	}

	x_I.resize(num_age);
	x_I[0]=N_pop[0]/birth_rate;
	for(int i=1; i<num_age; i++)
		x_I[i]=b[i-1]*N_pop[i]/N_pop[i-1];

	temp_pred.resize(num_age);

	double Sp=0.0;
	double am=age[num_age/2];
	for(int i=0; i<num_age; i++){
		Sp += N_pop[i];
		if(Sp>0.5){
			am=age[i+1];
			break;
		}
	}
	eta=v[1].size()==1 ? birth_rate : -1.0/am*log(1.0-Sp);
	
}

//***************************************************************
//***************************************************************