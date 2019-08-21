
//***************************************************************************
//***************************************************************************

class Event{
protected:
	double _exe_time;

public:
	Event() : _exe_time(-1E20){ }
	virtual void execute()=0;
	void schedule(Event_manager* event_manager, const double time_from_now=0);
	double exe_time() const{
		return _exe_time;
	}
	virtual bool event_empty() const{
		return false;
	}
	bool cancel(Event_manager* event_manager);
};

class Empty_event : public Event{
public:
	bool event_empty() const{
		return true;
	}
	void execute(){ }
};
//***************************************************************************
//***************************************************************************

class Move_state_event : public Event{
protected:
	Human* human;
public:
	Move_state_event(){}
	Move_state_event(Human* h) : human(h) {}
	void execute();
};

class Death_event : public Event{
protected:
	Human* human;
public:	
	Death_event(){}
	Death_event(Human* h) : human(h) {}
	void execute();
};

class Leave_itn_event : public Event{
protected:
	Human* human;
public:
	Leave_itn_event(){}
	Leave_itn_event(Human* h) : human(h) {}
	void execute();
};

class Pev_dose_event : public Event{
protected:
	Human* human;
public:
	Pev_dose_event(){}
	Pev_dose_event(Human* h) : human(h) {}
	void execute();
};

//***************************************************************************
//***************************************************************************

class Repeated_event : public Event{
protected:
	double interval;
	int max_rounds;
	int num_rounds;
	Village* village;
public:
	Repeated_event(Village* v, const double freq, const int m=10000, const double t=0);
	virtual void execute();
	virtual void event_occurs()=0;
	void skip_to(const double t);
};

//***************************************************************************
//***************************************************************************

class Output_event : public Repeated_event{
public:
	Output_event(Village* v, const double d, const double time_from_now=0) : Repeated_event(v, d, 1000000000, time_from_now) { }
	void event_occurs();
};

class Mosq_update_event : public Repeated_event{
protected:
public:
	Mosq_update_event(Village* v, const double d=0.25) : Repeated_event(v, d, 1000000000, 1E-7) { }
	void event_occurs();
};

class FOIM_update_event : public Repeated_event{
protected:
public:
	FOIM_update_event(Village* v, const double d=10.0) : Repeated_event(v, d, 1000000000, 3E-7) { }
	void event_occurs();
};

class Mass_birth_event : public Repeated_event{
protected:
public:
	Mass_birth_event(Village* v, const double d=5.0) : Repeated_event(v, d, 1000000000, 2E-7) { }
	void event_occurs();
};

//***************************************************************************
//***************************************************************************

template<class F>
class Repeated_intervention_event : public Repeated_event{
protected:
double prob;
	int target;
	double age0;
	double age1;
	F f;
	int corr_index;
	bool to_update_foim;

public:
	
	Repeated_intervention_event(Village* v, F f0, const string &prefix, 
		const int c, const double t=0, const bool foim=false, const string &suffix="");
	virtual void event_occurs();
};

//***************************************************************************
//***************************************************************************
template<class F>
class Single_intervention_event : public Event{
protected:
	Village* village;
	double p;
	F f;
public:
	Single_intervention_event(Village* v, F f0, const double p0, const double t);
	virtual void execute();
};

//***************************************************************************
//***************************************************************************

class Epi_event : public Event{
protected:
	Human* human;
public:
	Epi_event() : num_epi(0){ }
	Epi_event(Human* h) : num_epi(0), human(h) { }
	void execute();
	int num_epi;
};

//***************************************************************************
//***************************************************************************