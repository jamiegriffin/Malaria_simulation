//***************************************************************************
//***************************************************************************

class Event_manager{
private:
	double t;
	Calendar_queue Q;
	double span;
	double width;
	Empty_event empty_event;
	void execute_event();
	inline bool keep_going();
	
public:	
	Event_manager() : t(0) { } 
	void clear();
	void setup(const double t0, const double span0=2*365, const double width0=0.01);
	double time_now() const;

	void run(const double run_time);
	void run_det(const double run_time, ODE &ode);
	void add(Event* e);
	bool cancel(Event *e);
	bool find(Event *e);
	void jump_by(const double t0);
};
//***************************************************************************
//***************************************************************************
