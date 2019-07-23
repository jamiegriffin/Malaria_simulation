//***************************************************************************
//***************************************************************************

void Event_manager::setup(const double t0, const double span0, const double width0){
	t=t0;
	span=span0;
	width=width0;
	Q.setup(t, span, width);
}
double Event_manager::time_now() const{
	return t;
}

void Event_manager::jump_by(const double t0){
	t += t0;
}
void Event_manager::execute_event(){
	Event* e = Q.next_event();
	if(e!=0){
		if(e->exe_time()<t){
			cout << t << '\t' << e->exe_time() << '\t' 
				<< typeid(*e).name() << endl;
			error_crit("e->exe_time()<t in Event_manager::execute_event");
		}
		t = e->exe_time();
		e->execute();
	}
}
bool Event_manager::keep_going(){
	return !Q.next_event_empty();
}
void Event_manager::run(const double run_time){
	empty_event.schedule(this, run_time);
	const double end_time=empty_event.exe_time();
	while(keep_going())
		execute_event();
	empty_event.cancel(this);
	t=end_time;
}

void Event_manager::run_det(const double run_time, ODE &ode){
	empty_event.schedule(this, run_time);
	const double end_time=empty_event.exe_time();
	do{
		const double next_time=min(end_time, Q.next_time());
		if(next_time-t > 0){
			ode.step(next_time-t);
			t = next_time;
		}
		if(keep_going())
			execute_event();
	}while(keep_going());

	if(end_time-t > 0)
		ode.step(end_time-t);

	empty_event.cancel(this);
	t = end_time;
}

void Event_manager::add(Event* e){
	if(e->exe_time()<t){
		cout << t << '\t' << e->exe_time() << '\t' << typeid(*e).name() << endl;
		error_crit("e->exe_time()<t in Event_manager::add");
	}
	Q.add(e);
}

bool Event_manager::cancel(Event* e){
	return Q.remove(e);
}
bool Event_manager::find(Event* e){
	return Q.find(e);
}
void Event_manager::clear(){
	Q.clear();
	t=0;
}

//***************************************************************************
//***************************************************************************
