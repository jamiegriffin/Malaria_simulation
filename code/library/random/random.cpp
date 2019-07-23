//***************************************************************************
//***************************************************************************

double uniform(Ran &r) {
	double u;
	do {
		u = (r.ran64() >> 11)*1.1102230246251565e-16;
	} while (u==0.0 || u==1.0);
	return u;
}

uint64_t ran64(const uint64_t n, Ran &r){
	uint64_t j=0xffffffffffffffffULL;
	const uint64_t k=j/n;
	while(j/n==k){
		j=r.ran64();
	}
	return j % n;
}
unsigned int ran_int(const unsigned int n, Ran &r){
	return static_cast<unsigned int>(ran64(n, r));
}

double rnorm(Ran &r) {
	return invnorm(uniform(r));
}

int rbinom2(const double p, const int n, Ran &r) {
	if (n>20  && n*p<=7.0) {
		// Devroye waiting time method 2
		int x = 0;
		const double q = -log(1.0-p);
		double S = 0.0;
		do {
			if (x==n)
				return n;
			S += -log(uniform(r))/(n-x);
			x++;
		} while (S<=q);
		return x-1;
	}
	else if (n<=80) {
		int x = 0;
		const uint64_t m = static_cast<uint64_t>(p*double(0xffffffffffffffffULL));
		for (int i = 0; i<n; i++)
			if (r.ran64() < m)
				x++;
		return x;
	}
	else {
		// Devroye recursive binomial generator
		const double a = floor((n+1.0)*p);
		const double y = rbeta(a, n+1.0-a, r);
		return y<=p ? a + rbinom((p-y)/(1.0-y), n-a, r)	: rbinom(p/y, a-1, r);
	}
}
int rbinom(const double p, const int n, Ran &r) {
	if (n==0 || p==0)
		return 0;
	else if (p==1)
		return n;
	const double p1 = p<=0.5 ? p : 1.0-p;
	const int x = rbinom2(p1, n, r);
	return p1==p ? x : n-x;
}

int rpoisson(const double mu0, Ran &r) {
	double mu = mu0;
	int x = 0;
	while (mu > 10.0) {
		const int n = static_cast<int>(floor(mu*0.875));
		const double g = rgamma(n, r);
		if (g > mu) {
			return x + rbinom(mu/g, n-1, r);
		}
		else {
			x += n;
			mu -= g;
		}
	}
                                                        
	const double q = exp(-mu);
	double y = 1.0;
	do {
		y *= uniform(r);
		x++;
	} while (y > q);
	return x-1;
}

void rmnorm(const vector<vector<double> > &C, vector<double> &z, Ran &r){
	const int n=z.size();
	for(int i=0; i<n; i++)
		z[i]=0.0;
	for(int j=0; j<n; j++){
		const double x=rnorm(r);
		for(int i=j; i<n; i++)
			z[i] += C[i][j]*x;
	}
}

double rgamma(const double a, Ran &r){
	if(a==1.0)
		return -log(uniform(r));
	else if(a==2.0)
		return -log(uniform(r)*uniform(r));
	else if(a==3.0)
		return -log(uniform(r)*uniform(r)*uniform(r));
	else if(a==4.0)
		return -log(uniform(r)*uniform(r)*uniform(r)*uniform(r));

	else if(a>1.0){
		// Marsaglia and Tsang 2000
		const double d=a-0.333333333333333333;
		const double b=-sqrt(9.0*d);
		const double c = -1.0/b;
		double v;
		while(true){
			double x;
			do{
				x=rnorm(r);
			}while (x <= b);
			v = 1.0 + c*x;
			v *= v*v;
			const double u = uniform(r);
			const double y = x*x;
			if (u < 1.0 - 0.331*y*y)
				break;
			if (log(u) < 0.5*y + d*(1.0-v+log(v)))
				break;
		}
		return d*v;
	}
	else if (a<1.0 && a>0.0){
		// Non-Uniform Random Variate Generation Luc Devroye, p304
		const double a1=1.0/a;
		while(true){
			const double u=uniform(r);
			const double e0=-log(uniform(r));
			if(u<1.0-a){
				const double x=pow(u, a1);
				if(x<e0)
					return x;
			}
			else{
				const double e=-log((1.0-u)/a);
				const double x=pow(1.0-a+a*e, a1);
				if(x<e0+e)
					return x;
			}
		}
	}
	else{
		error_crit("a must be >0 in rgamma()");
		return 0.0;
	}
}

double rbeta(const double a, const double b, Ran &r) {
	const double x = rgamma(a, r);
	const double y = rgamma(b, r);
	return x/(x+y);
}

//***************************************************************************
//***************************************************************************

int discrete_dev(const vector<double> &p, const double u0, const bool norm){
	double S=0.0;
	const double u=(norm ? u0 : u0*sum(p));
	const int n=p.size();
	for(int j=0; j<n; j++){
		if(p[j]<-1.0E-20)
			error_crit("All elements of input vector in discrete_dev() must be non-negative");
		S += p[j];
		if(u<S) return j;
	}
	return n-1;
}

int discrete_dev(const vector<double> &p, Ran &r, const bool norm){
  	return discrete_dev(p, uniform(r), norm);
}

void Guide_table::setup(const vector<double> &p){
	for(int i=0; i<p.size(); i++)
		if(p[i]<-1.0E-12)
			error_crit("All elements of input vector in Guide_table must be non-negative");
	M=p.size()*a;
	if(q.size()!=p.size()) q.resize(p.size());
	g.assign(M, 0);
	const double S=sum(p);
	for(int i=0; i<p.size()-1; i++){
		q[i]=(i==0 ? 0.0 : q[i-1]) + p[i]/S;
		g[min(int(floor(q[i]*M)), M-1)]=i+1;
	}
	q.back()=1.0;
	for(int j=1; j<M; j++)
		g[j]=max(g[j-1], g[j]);
}
Guide_table::Guide_table(const vector<double> &p, const int aa) : a(aa){
	setup(p);
}
int Guide_table::operator()(Ran &r) const{
	return (*this)(uniform(r));
}
int Guide_table::operator()(const double u) const{
	if(M==0) error_crit("guide table not yet set up");
	if(u<0.0 || u>1.0) error_crit("u out of range in Guide_table::operator())");
	int i=g[min(int(u*M), M-1)];
	while(i>0 && q[i-1]>u) i--;
	return i;
}

//***************************************************************************
//***************************************************************************
