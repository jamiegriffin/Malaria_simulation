//***************************************************************************
//***************************************************************************

const int num_mosq=3;
const double dy=365;

//***************************************************************************
//***************************************************************************

void save_state(const vector<double> &v, ofstream &out){
	const int n=v.size();
	output_binary(n, out);
	for(int i=0; i<n; i++)
		output_binary(v[i], out);
}
void use_state(vector<double> &v, ifstream &in){
	int n;
	input_binary(n, in);
	v.resize(n);
	for(int i=0; i<n; i++)
		input_binary(v[i], in);
}

struct Output{
	bool newline;
	bool binary;
	ofstream& out;
	Output(ofstream &out0, const bool b) : out(out0), binary(b), newline(true){	}
	void operator()(const double x){
		if(binary)
			output_binary(x, out);
		else{
			if(!newline)
				out << '\t';
			out << x;
			newline=false;
		}
	}
	void operator()(const string &s){
		if(!binary){
			if(!newline)
				out << '\t';
			out << s;
			newline=false;
		}
	}
	void endl(){
		if(!binary){
			out << '\n';
			newline=true;
		}
	}
};

struct Output_var{
	enum Quantity {prev, pcr_prev, clin_inc, sev_inc, EIR, prop, inf_inc, itn_cov, SM};
	Quantity quantity;
	double age0;
	double age1;
	Output_var() {}
	Output_var(const string &q, const double a0, const double a1){
		age0=a0>=0 ? a0*dy : 0;
		age1=a1>=0 ? a1*dy : 200*dy;

		if(q=="prev")
			quantity=prev;
		else if(q=="pcr_prev")
			quantity=pcr_prev;
		else if(q=="clin_inc")
			quantity=clin_inc;
		else if(q=="inf_inc")
			quantity=inf_inc;
		else if(q=="sev_inc")
			quantity=sev_inc;
		else if(q=="EIR")
			quantity=EIR;
		else if(q=="prop")
			quantity=prop;
		else if(q=="itn_cov")
			quantity=itn_cov;
		else if (q=="SM")
			quantity = SM;
		else
			error_crit("Output quantity " + q + " not implemented");
	}
};

class Event_manager;
class Event;

struct Simulation;

struct Village;
struct Village_het;
struct Human;
struct ITN_IRS;
struct Mosq_pop;
struct Parms;
struct Human_parms;

struct Drug{
	double dur_T;	
	double dur_P;
	double shape_P;
	double rc;		
	double efficacy;

	Guide_table age_guide;
	double max_age;
	double max_day;
	double dt;
	bool p_protect;
	
	vector<vector<double> > p_infect;
	
	double prob_infected(const double time_since_trt, const double age){
		if(!p_protect)
			return 1.0;
		if(time_since_trt < 0)
			error_crit("time_since_trt < 0, " + to_string(time_since_trt));
		return (time_since_trt>max_day-1E-12) ? 1.0 : p_infect[find_age_group(age)][time_since_trt/dt];
	}

	int find_age_group(const double age1) const{
		return !p_protect ? 0 : ((age1>=max_age ? p_infect.size()-1 : age_guide(age1/max_age)));
	}
	double choose_dur(Ran& ran) const{
		if(p_protect)
			error_crit("Drug::choose_dur() should only be called when p_protect is false");
		return dur_P/shape_P*rgamma(shape_P, ran);
	}
	Drug() : dur_T(5), dur_P(25), shape_P(1.0), rc(0.32), efficacy(0.75) { } 
	
	void setup(const double dT, const double r, const double e, const string &file){
	
		p_protect=from_map_bool("p_protect");

		dur_T=dT;
		rc=r;
		efficacy=e;

		if(p_protect){
			vector<vector<double> > v;
			import_file(file, v, false);
			transpose(v);
			vector<double> ages=v[1];
			for(int a=0; a<ages.size(); a++)
				ages[a] *= 365.0;
			sort(ages.begin(), ages.end());
			vector<double>::iterator it=unique(ages.begin(), ages.end());
			ages.erase(it, ages.end());
			max_age=max_el(v[2])*365.0;
			const int na=ages.size();
			for(int a=0; a<na; a++)
				ages[a]=(a<na-1 ? ages[a+1] : max_age) - ages[a];

			age_guide.setup(ages);
			
			max_day=max_el(v[0]);
			dt=v[0][1]-v[0][0];
			const int nd=v[0].size()/na;
			p_infect.assign(na, vector<double>(nd));
			for(int a=0; a<na; a++)
				for(int d=0; d<nd; d++)
					p_infect[a][d]=1.0-(d==nd-1 ? v[3][a*nd + d] : 0.5*(v[3][a*nd + d]+v[3][a*nd + d + 1]));

			double Sp=0;
			for(int d=0; d<nd; d++)
				for(int a=0; a<na; a++)
					Sp += (1.0-p_infect[a][d])*dt;
			Sp /= na;
			dur_P=Sp;
		}
	} 
};

struct ITN_IRS{
	double itn_time;
	double irs_time;
	bool itn;
	bool irs;
	ITN_IRS() : itn(false), irs(false), itn_time(1E20), irs_time(1E20) { }
	ITN_IRS(const bool i1, const bool i2, const double t1, const double t2) : itn(i1), irs(i2),itn_time(t1), irs_time(t2) { }

	void save_state(ofstream &out) {
		output_binary(itn_time, out);
		output_binary(irs_time, out);
		output_binary(itn, out);
		output_binary(irs, out);
	}
	void use_state(ifstream &in) {
		input_binary(itn_time, in);
		input_binary(irs_time, in);
		input_binary(itn, in);
		input_binary(irs, in);
	}
};

#pragma warning(disable : 4996)
string extract_time_info(const char* format){
	char c[40];
	time_t tim=time(0);
	strftime(c, 40, format, localtime(&tim));
	return string(c);
}
string timestr(){
	return extract_time_info("%H.%M.%S");
}
string datestr(){
	return extract_time_info("%A %d %b %Y");
}

#include "equilibrium/equilibrium.h"
#include "equilibrium/util.h"

#include "equilibrium/malaria state.h"
#include "equilibrium/ode.h"
#include "equilibrium/malaria_parms.h"

#include "equilibrium/ode.cpp"

#include "events/event.h"
#include "events/queues.h"
#include "events/event_manager.h"

//***************************************************************************
//***************************************************************************