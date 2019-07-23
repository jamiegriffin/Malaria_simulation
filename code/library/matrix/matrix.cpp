
void cholesky(const vector<vector<double> > &A, vector<vector<double> > &C, const double eps){
	C = A;
	const int n=A.size();
	if(eps>0.0)
		for(int i=0; i<n; i++)
			C[i][i] *= 1.0+eps;
	
	for (int i=0;i<n;i++) {
		for (int j=i;j<n;j++) {
			double S=C[i][j];
			for (int k=i-1;k>=0;k--)
				S -= C[i][k]*C[j][k];
			C[j][i]= i==j ? sqrt(S) : S/C[i][i];
		}
	}	
	for (int i=0;i<n;i++) {
		for (int j=i+1;j<n;j++) {
			C[i][j]=0.0;
		}
	}
}
