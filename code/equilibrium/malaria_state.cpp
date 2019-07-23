
//***************************************************************
//***************************************************************

void Malaria_state::delete_storage(){
	if(num_het>=0){
		delete [] y;
		y=0;
		num_het=-1;
		num_age=-1;
		human_state.clear();
		mosq_state.clear();
	}
}


void Malaria_state::setup(const Malaria_state &z){
	if(num_age != z.num_age || num_het != z.num_het)
		setup(z.num_age, z.num_het);
	copy(z.y, z.y+total_size, y);
}
void Malaria_state::setup(const int na, const int nh){

	delete_storage();
	num_het=nh;
	num_age=na;
	human_state.assign(2, vector<Human_state>(num_het, num_age));		

	const int mosq0_size=Mosq_state::size();
	const int mosq_size=num_mosq*mosq0_size;
	const int human0_size=human_state[0][0].size();
	const int human_size=num_het*human0_size;
	total_size=mosq_size + 2*human_size;
	init_size=mosq_size + human_size;
	current_size=init_size;
	
	y=new double[total_size];

	for(int m=0; m<num_mosq; m++)
		mosq_state.push_back(Mosq_state(y + m*mosq0_size));
	for(int k=0; k<2; k++)
		for(int j=0; j<num_het; j++)
			human_state[k][j].setup(y + mosq_size + (num_het*k + j)*human0_size);
}

void Malaria_state::setup(const Parms &parms){

	if(num_age != parms.num_age || num_het != parms.num_het)
		setup(parms.num_age, parms.num_het);
	else{
		for(int l=0; l<total_size; l++)
			y[l]=0;
		for(int m=0; m<num_mosq; m++)
			mosq_state[m].set_init();
	}

	for(int j=0; j<num_het; j++)
		for(int i=0; i<num_age; i++)
			human_state[0][j].S[i]=parms.het_wt[j]*parms.N_pop[i];
	
}

void Malaria_state::add_intervention(const double coverage){
	current_size=total_size;
	for(int j=0; j<num_het; j++){
		for(int l=0; l<human_state[0][j].num_prev_states(); l++){
			for(int i=0; i<num_age; i++){
				const double p=human_state[0][j].state(l)[i];
				human_state[0][j].state(l)[i]=(1-coverage)*p;
				human_state[1][j].state(l)[i]=coverage*p;
			}
		}
		for(int l=human_state[0][j].num_prev_states(); l<human_state[0][j].num_states(); l++)
			for(int i=0; i<num_age; i++)
				human_state[1][j].state(l)[i]=human_state[0][j].state(l)[i];
	}
}

void Malaria_state::operator = (const Malaria_state &z){
	setup(z);
}

//***************************************************************
//***************************************************************