/*********************************************************************
*********************************************************************/

double sum(const double* y, const int n){
	double sum=0.0;
	for(int i=0; i<n; i++){
		sum += y[i];
	}
	return sum;
}

double mean(const double* y, const int n){
	return sum(y, n)/double(n);
}
double sum(const vector<double> &v){
	double sum=0.0;
	const int n=v.size();
	for(int i=0; i<n; i++){
		sum += v[i];
	}
	return sum;
}

double mean(const vector<double> &v){
	return sum(v)/static_cast<double>(v.size());
}

/*********************************************************************
*********************************************************************/
