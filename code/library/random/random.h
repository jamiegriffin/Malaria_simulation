//*********************************************************************
//*********************************************************************

class Xor64star {
public:
	Xor64star() : x(0ULL) { }
	Xor64star(const uint64_t seed){
		setup(seed);
	}
	void setup(const uint64_t seed) {
		x = seed;
		do {
			x = x*3935559000370003845ULL + 2691343689449507681ULL;
		} while (x==0ULL);
		ran64();
	}

	uint64_t ran64() {
		// Constants from Vigna 2016
		// Method based on Marsaglia 2003
		x ^= (x >> 12);
		x ^= (x << 25);
		x ^= (x >> 27);
		return x*2685821657736338717ULL;
	}
	void save_state(ofstream &out) {
		output_binary(x, out);
	}
	void use_state(ifstream &in) {
		input_binary(x, in);
	}

private:
	uint64_t x;
};

typedef Xor64star Ran;

//*********************************************************************
//*********************************************************************

double uniform(Ran &r);
uint64_t ran64(const uint64_t n, Ran &r);
unsigned int ran_int(const unsigned int n, Ran &r);
double rnorm(Ran &r);
int rbinom(const double p0, const int n, Ran &r);
int rpoisson(const double mu0, Ran &r);
double rgamma(const double a, Ran& r);
double rbeta(const double a, const double b, Ran& r);

void rmnorm(const vector<vector<double> > &C, vector<double> &z, Ran& r);

int discrete_dev(const vector<double> &p, const double u0, const bool norm = false);
int discrete_dev(const vector<double> &p, Ran& r, const bool norm = false);

class Guide_table {
private:
	vector<double> q;
	vector<int> g;
	int a;
	int M;
public:
	Guide_table(const int aa = 5) : a(aa), M(0) { };
	Guide_table(const vector<double> &, const int aa = 5);
	void setup(const vector<double> &p);
	int operator()(Ran &r) const;
	int operator()(const double u) const;
};

//*********************************************************************
//*********************************************************************





