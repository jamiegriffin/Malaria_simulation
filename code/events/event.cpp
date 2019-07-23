//***************************************************************************
//***************************************************************************

void Event::schedule(Event_manager* event_manager, const double time_from_now){
	if(time_from_now < 0){
		cout << event_manager->time_now() << '\t' <<  time_from_now << '\t' << typeid(*this).name() << '\n';
		error_crit("time_from_now < 0 in Event::schedule");
	}
	_exe_time=time_from_now + event_manager->time_now();
	event_manager->add(this);
}
bool Event::cancel(Event_manager* event_manager){
	return event_manager->cancel(this);
}

//***************************************************************************
//***************************************************************************

void Move_state_event::execute(){
	human->move_state();
}
void Death_event::execute(){
	human->die();
}
void Leave_itn_event::execute(){
	human->leave_itn();
}

//***************************************************************************
//***************************************************************************

Repeated_event::Repeated_event(Village* v, const double freq, const int m, const double t) : village(v), interval(freq), num_rounds(0), max_rounds(m) {
	schedule(&village->simulation->event_manager, t);
}
void Repeated_event::execute(){
	event_occurs();
	if(++num_rounds < max_rounds)
		schedule(&village->simulation->event_manager, interval);
}
void Repeated_event::skip_to(const double t){
	double t1=_exe_time;
	while(t1<t){
		t1 += interval;
		num_rounds++;	
	}
	if(num_rounds<max_rounds){
		const double t0=village->simulation->event_manager.time_now();
		schedule(&village->simulation->event_manager, t1-t0);
	}
}

//***************************************************************************
//***************************************************************************

void Output_event::event_occurs(){
	village->output(num_rounds-1, _exe_time, interval);
}

void Mosq_update_event::event_occurs(){
	village->update_mosq(interval);
}

void FOIM_update_event::event_occurs(){
	village->update_foim();
}

void Mass_birth_event::event_occurs(){
	village->mass_birth();
}

//***************************************************************************
//***************************************************************************

template<class F>
Repeated_intervention_event<F>::Repeated_intervention_event(Village* v, F f0, const string &prefix, 
		const int c, const double t, const bool foim, const string &suffix) 
		: Repeated_event(v,
				(prefix=="itn" ? dy : from_map_suffix(prefix + "_frequency", suffix, 0.0, 10000.0, 1.0)*dy),
				from_map_int_suffix(prefix + "_max_rounds", suffix, 0, 2000000000, 10000),
				t + from_map_suffix(prefix + "_start", suffix, 0.0, 10000.0, 0.0)*dy), 
		f(f0), corr_index(c), to_update_foim(foim) {
	prob=from_map_suffix(prefix + "_coverage", suffix, 0, 1);
	age0=from_map_suffix(prefix + "_age0", suffix, 0.0, 200.0, 0.0)*dy;
	age1=from_map_suffix(prefix + "_age1", suffix, 0.0, 200.0, 199.0)*dy;		
	target=from_map_int_suffix(prefix + "_target", suffix, 0, 3, 0);
	if(target>1 && v->village_het.size()<10)
		error_crit("Need at least 10 heterogeneity groups (num_het) if target = 2 or 3");

}
template<class F>
void Repeated_intervention_event<F>::event_occurs(){
	if(!is_same<F, Itn_pulse_func>::value || village->simulation->itn_now())
		village->pulse(prob, target, corr_index, f, age0, age1, to_update_foim);
}

template <class F>
void repeated_intervention_event(Village* v, F f0, const string &prefix, 
		const int c, const double t=0, const bool foim=false, const string &suffix=""){
	v->simulation->events.push_back(new Repeated_intervention_event<F>(v, f0, prefix, c, t, foim, suffix));
}

//***************************************************************************
//***************************************************************************
template <class F>
Single_intervention_event<F>::Single_intervention_event(Village* v, F f0, const double p0, const double t) : village(v), f(f0), p(p0) {
	schedule(&village->simulation->event_manager, t);
}
template <class F>
void Single_intervention_event<F>::execute(){
	f(village, p);
}
template <class F>
void single_intervention_event(Village* v, F f0, const double p0, const double t){
	v->simulation->events.push_back(new Single_intervention_event<F>(v, f0, p0, t));
}

//***************************************************************************
//***************************************************************************

void Epi_event::execute(){
	human->give_epi(num_epi);
	num_epi++;
	vector<double>* epi_ages=&human->village_het->village->simulation->epi_ages;
	if(num_epi<epi_ages->size())
		schedule(&human->village_het->village->simulation->event_manager, (*epi_ages)[num_epi]-(*epi_ages)[num_epi-1]);
}

//***************************************************************************
//***************************************************************************