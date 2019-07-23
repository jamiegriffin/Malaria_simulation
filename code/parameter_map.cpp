//***************************************************************************
//***************************************************************************

#include <map>

typedef map<string, pair<double, string> > Parameter_map;
Parameter_map parameter_map;
void add_to_map(const string &s, const double x){
	parameter_map[s]=pair<double, string>(x, "");
}
void modify_map(const string &s, const double x){
	Parameter_map::iterator p=parameter_map.find(s);
	if(p == parameter_map.end())
		error_crit("Trying to modify parameter "+s+" from the command line, but that parameter name was not found in any of the input files");
	p->second.first=x;
}
bool in_map(const string &s){
	return parameter_map.find(s) != parameter_map.end();
}

void import_to_map(const string &file, const bool add=false){
	ifstream input(file.c_str(), ios::in);
	if(!input)
		error_crit("Can't open file "+file+" for input");
	while(!input.eof()){
		string s;
		double x;
		input >> s >> x;
		if(!input.eof()){
			if(!input.good())
				error_crit("Error in importing data from "+file);
			if(add){
				if(in_map(s))
					modify_map(s, x);
				else
					parameter_map.insert(pair<string, pair<double, string> >(s, pair<double, string>(x, file)));
			}
			else{
				pair<Parameter_map::const_iterator, bool> pb=parameter_map.insert(pair<string, pair<double, string> >(s, pair<double, string>(x, file)));
				if(!pb.second)
					error_crit("Importing parameter "+s+" from file\n"+file+", but a parameter with that name has already been imported from file\n"+pb.first->second.second);	
			}
		}
	}
	input.close();
}

double from_map(const string &s, const double low=0, const double high=1E20, const double x=-99999){
	Parameter_map::const_iterator p=parameter_map.find(s);
	const bool found=p!=parameter_map.end();
	if(!found && x==-99999)
		error_crit("Parameter "+s+" was not found in any of the input files");

	const double y=found ? p->second.first : x;
	const double eps=1E-8;
	if(y<low && y>low-eps)
		return low;
	else if(y>high && y<high+eps)
		return high;
	else if(y<=low-eps || y>=high+eps)
		error_crit("Parameter "+s+" has value "+to_string(y)+" , should be between "+to_string(low)+" and "+to_string(high)+ (found ? " in file \n"+p->second.second : " as default"));

	return y;
}
long long from_map_ll(const string &s, const double low=0, const double high=9223372036854775807.0, const double x=-99999){
	Parameter_map::const_iterator p=parameter_map.find(s);
	const bool found=p!=parameter_map.end();
	if(!found && x==-99999)
		error_crit("Parameter "+s+" was not found in any of the input files");

	const double y=found ? p->second.first : x;

	if(y != floor(y+0.5))
		error_crit("Parameter "+s+" has value "+to_string(y)+" , should be an integer" + (found ? " in file \n"+p->second.second : " as default"));
	if(y<low || y>high)
		error_crit("Parameter "+s+" has value "+to_string(y)+" , should be between "+to_string(low)+" and "+to_string(high)+ (found ? " in file \n"+p->second.second : " as default"));
	return y;
}
int from_map_int(const string &s, const double low=0, const double high=2147483647.0, const double x=-99999){
	return from_map_ll(s, low, high, x);
}
bool from_map_bool(const string &s, const double x=-99999){
	Parameter_map::const_iterator p=parameter_map.find(s);
	const bool found=p!=parameter_map.end();
	if(!found && x==-99999)
		error_crit("Parameter "+s+" was not found in any of the input files");

	const double y=found ? p->second.first : x;

	if(y==0)
		return false;
	else if(y==1)
		return true;
	else
		error_crit("Parameter "+s+" has value "+to_string(y)+" , only 0 or 1 allowed" + (found ? " in file \n"+p->second.second : " as default"));

	return true;
}

double from_map_suffix(const string &s, const string &suffix, const double low=0, const double high=1E20, const double x=-99999){
	return in_map(s + suffix) ? from_map(s + suffix, low, high, x) : from_map(s, low, high, x);
}
long long from_map_ll_suffix(const string &s, const string &suffix, const double low=0, const double high=9223372036854775807.0, const double x=-99999){
	return in_map(s + suffix) ? from_map_ll(s + suffix, low, high, x) : from_map_ll(s, low, high, x);
}
int from_map_int_suffix(const string &s, const string &suffix, const double low=0, const double high=2147483647.0, const double x=-99999){
	return from_map_ll_suffix(s, suffix, low, high, x);
}
bool from_map_bool_suffix(const string &s, const string &suffix, const double x=-99999){
	return in_map(s + suffix) ? from_map_bool(s + suffix, x) : from_map_bool(s, x);
}

//***************************************************************************
//***************************************************************************