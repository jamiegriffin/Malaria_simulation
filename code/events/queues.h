//***************************************************************************
//***************************************************************************

template<class T>
void shrink_vector(vector<T> &v){
	vector<T>(v).swap(v);
}

struct later_event : public binary_function<Event, Event, bool>{
	bool operator()(const Event *const e1, const Event *const e2) const{
		if(e1->exe_time()==e2->exe_time()){
			if(e1->event_empty() && !e2->event_empty())
				return true;
			else if(!e1->event_empty() && e2->event_empty())
				return false;
			else
				return false;
		}
		else
			return (e1->exe_time() > e2->exe_time());
	}
};

//***************************************************************************
//***************************************************************************
class Calendar_queue{
	// A calendar queue implemented as a vector of vectors, 
	// with the current bucket from which events are being executed 
	// kept as a priority queue, and the other buckets completely unsorted

private:
	long long current_bucket;
	long long current_virtual_bucket;
	long long N_events;
	double width;
	vector<vector<Event*> > Q;
	long long nq;
	
	long long find_virtual_bucket(Event* e) const;
	void update_current_bucket();

public:
	Calendar_queue() 
			: current_bucket(0), current_virtual_bucket(0), N_events(0) {}

	void setup(const double t, const double span=2*365, const double w=0.1){
		width=w;
		nq=static_cast<long long>(span/w+0.5);
		Q.resize(nq);
		current_virtual_bucket=static_cast<long long>(floor(t/width));
	}
	bool empty() const{
		if(N_events==0)
			cout << "N_events==0\n";
		return (N_events==0);
	}
	long long size() const{
		return N_events;
	}
	void add(Event* e);
	Event* next_event();
	double next_time();
	bool next_event_empty();
	bool remove(Event* e);
	bool find(Event* e);
	void clear();
};

long long Calendar_queue::find_virtual_bucket(Event* e) const{
	return static_cast<long long>(floor(e->exe_time()/width));
}
void Calendar_queue::update_current_bucket(){
	bool advance_bucket=Q[current_bucket].empty() || 
			(find_virtual_bucket(Q[current_bucket].front()) > current_virtual_bucket);
	//if(advance_bucket)
	//	shrink_vector(Q[current_bucket]); // shrink bucket to fit its remaining contents before advancing to next one - saves memory, but may be slower
	while(advance_bucket){
		current_virtual_bucket++;
		current_bucket= modulo(current_virtual_bucket, nq);
		if(!Q[current_bucket].empty()){
			// make the new current bucket into a priority queue
			make_heap(Q[current_bucket].begin(), Q[current_bucket].end(), later_event());
			advance_bucket= (find_virtual_bucket(Q[current_bucket].front()) > current_virtual_bucket);
		}
	}
}
void Calendar_queue::add(Event* e){ 
	const long long vb=find_virtual_bucket(e);
	if(vb<current_virtual_bucket)
		error_crit("virtual_bucket<current_virtual_bucket in Calendar_queue::add");
	const long long bucket=modulo(vb, nq);
	Q[bucket].push_back(e);
	// if the bucket is the current bucket, use push_heap to maintain the priority queue ordering
	if(bucket == current_bucket) 
		push_heap(Q[bucket].begin(), Q[bucket].end(), later_event());
	N_events++;
}

Event* Calendar_queue::next_event(){ 
	if(N_events==0LL)
		return 0;
	update_current_bucket();
	Event* e=Q[current_bucket].front();
	pop_heap(Q[current_bucket].begin(), Q[current_bucket].end(), later_event());
	Q[current_bucket].pop_back();
	N_events--;
	return e;
}
double Calendar_queue::next_time(){
	if(N_events==0LL)
		return 1E20;
	update_current_bucket();
	return Q[current_bucket].front()->exe_time();
}
bool Calendar_queue::next_event_empty(){
	if(N_events==0LL)
		return true;
	update_current_bucket();
	return Q[current_bucket].front()->event_empty();
}

void Calendar_queue::clear(){
	for(int i=0; i<Q.size(); i++)
		Q[i].clear();
	current_bucket=0LL;
	current_virtual_bucket=0LL;
	N_events=0LL;
}

bool Calendar_queue::find(Event* e){
	const long long bucket=modulo(find_virtual_bucket(e), nq);
	return (::find(Q[bucket].begin(), Q[bucket].end(), e) != Q[bucket].end());
}
bool Calendar_queue::remove(Event* e){
	const long long bucket=modulo(find_virtual_bucket(e), nq);
	const bool removed=remove_once(Q[bucket], e);
	if(removed){
		if(bucket==current_bucket && Q[bucket].size()>0)
			make_heap(Q[bucket].begin(), Q[bucket].end(), later_event());
		N_events--;
		return true;
	}
	return false;
}

//***************************************************************************
//***************************************************************************

