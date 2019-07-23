//*********************************************************************
//*********************************************************************

void error_crit(const string &s){
	cerr << endl << s << endl;
	exit(1);
}

template<class T>
void output_binary(const T x, ofstream &out) {
	out.write((char *)&x, sizeof(T));
}
template<class T>
void input_binary(T &x, ifstream &in) {
	in.read((char *)&x, sizeof(T));
}

template<class T>
T modulo(const T x, const T y){
	const T r = x % y;
    return r<0 ? r+y : r;
}
double invlogit(const double x){
	const double w=exp(x);
	return w/(1+w);
}

//*********************************************************************
//*********************************************************************
