//***************************************************************************
//***************************************************************************

size_t count_vars(ifstream &input){
	size_t num_spaces=0;
	size_t num_c=0;
	bool end_line=false;
	while(!end_line){
		char c;
		input.get(c);
		num_c++;
		if(c == '\t' || c  == ' ')
			num_spaces++;
		end_line=(c == '\n');
	}
	for(size_t k=0; k<num_c; k++)
		input.unget();
	return num_spaces+1;
}

//***************************************************************************
//***************************************************************************
template<class T>
void transpose(vector<vector<T> > &v){
	vector<vector<T> > w;
	w.assign(v[0].size(), vector<T>(v.size()));
	for(size_t i=0; i<v.size(); i++)
		for(size_t j=0; j<v[0].size(); j++)
			w[j][i]=v[i][j];
	v=w;
}

template<class T>
void import_file(const string &file, vector<vector<T> > &v, const bool header=false){
	ifstream s(file.c_str(), ios::in);
	if(!s)
		error_crit("Can't open file "+file+" for input");
	const size_t n=count_vars(s);
	
	vector<T> x(n);
	if(header){
		string name;
		for(size_t i=0; i<n; i++)
			s >> name;
	}
	while(!s.eof()){
		for(size_t i=0; i<n; i++)
			s >> x[i];
		if(!s.eof()){
			if(!s.good())
				error_crit("Error in importing data from "+file);
			v.push_back(x);
		}
	}
	s.close();
}
//***************************************************************************
//***************************************************************************
template<class T1>
void import_file(const string &file, vector<T1> &v1){
	if(v1.size() != 0)
		error_crit("vector should be empty in import_file, importing from file "+file); 
	ifstream s(file.c_str(), ios::in);
	if(!s)
		error_crit("Can't open file "+file+" for input");
	T1 x1;
	while(!s.eof()){
		s >> x1;
		if(!s.eof()){
			if(!s.good())
				error_crit("Error in importing data from "+file);
			v1.push_back(x1);
		}
	}
	s.close();
}
template<class T1, class T2>
void import_file(const string &file, vector<T1> &v1, vector<T2> &v2){
	if(v1.size() != 0 || v2.size() != 0)
		error_crit("vector should be empty in import_file, importing from file "+file);
	ifstream s(file.c_str(), ios::in);
	if(!s)
		error_crit("Can't open file "+file+" for input");
	T1 x1;
	T2 x2;
	while(!s.eof()){
		s >> x1 >> x2;
		if(!s.eof()){
			if(!s.good())
				error_crit("Error in importing data from "+file);
			v1.push_back(x1);
			v2.push_back(x2);
		}
	}
	s.close();
}
template<class T1, class T2, class T3>
void import_file(const string &file, vector<T1> &v1, vector<T2> &v2, vector<T3> &v3){
	if(v1.size() != 0 || v2.size() != 0 || v3.size() != 0)
		error_crit("vector should be empty in import_file, importing from file "+file);
	ifstream s(file.c_str(), ios::in);
	if(!s)
		error_crit("Can't open file "+file+" for input");
	T1 x1;
	T2 x2;
	T3 x3;
	while(!s.eof()){
		s >> x1 >> x2 >> x3;
		if(!s.eof()){
			if(!s.good())
				error_crit("Error in importing data from "+file);
			v1.push_back(x1);
			v2.push_back(x2);
			v3.push_back(x3);
		}
	}
	s.close();
}
template<class T1, class T2, class T3, class T4>
void import_file(const string &file, vector<T1> &v1, vector<T2> &v2, vector<T3> &v3, vector<T4> &v4){
	if(v1.size() != 0 || v2.size() != 0 || v3.size() != 0 || v4.size() != 0)
		error_crit("vector should be empty in import_file, importing from file "+file);
	ifstream s(file.c_str(), ios::in);
	if(!s)
		error_crit("Can't open file "+file+" for input");
	T1 x1;
	T2 x2;
	T3 x3;
	T4 x4;
	while(!s.eof()){
		s >> x1 >> x2 >> x3 >> x4;
		if(!s.eof()){
			if(!s.good())
				error_crit("Error in importing data from "+file);
			v1.push_back(x1);
			v2.push_back(x2);
			v3.push_back(x3);
			v4.push_back(x4);
		}
	}
	s.close();
}
template<class T1, class T2, class T3, class T4, class T5>
void import_file(const string &file, vector<T1> &v1, vector<T2> &v2, vector<T3> &v3, vector<T4> &v4, vector<T5> &v5){
	if(v1.size() != 0 || v2.size() != 0 || v3.size() != 0 || v4.size() != 0 || v5.size() != 0)
		error_crit("vector should be empty in import_file, importing from file "+file);
	ifstream s(file.c_str(), ios::in);
	if(!s)
		error_crit("Can't open file "+file+" for input");
	T1 x1;
	T2 x2;
	T3 x3;
	T4 x4;
	T5 x5;
	while(!s.eof()){
		s >> x1 >> x2 >> x3 >> x4 >> x5;
		if(!s.eof()){
			if(!s.good())
				error_crit("Error in importing data from "+file);
			v1.push_back(x1);
			v2.push_back(x2);
			v3.push_back(x3);
			v4.push_back(x4);
			v5.push_back(x5);
		}
	}
	s.close();
}
//***************************************************************************
//***************************************************************************
template<class T1>
void export_file(const string &file, const vector<T1> &v1){
	ofstream s(file.c_str(), ios::out);
	if(!s)
		error_crit("Can't create file "+file+" for output");
	for(size_t i=0; i<v1.size(); i++)
		s << v1[i] << '\n';
	s.close();
}
template<class T1, class T2>
void export_file(const string &file, const vector<T1> &v1, const vector<T2> &v2){
	ofstream s(file.c_str(), ios::out);
	if(!s)
		error_crit("Can't create file "+file+" for output");
	if(v2.size() != v1.size())
		error_crit("vectors should be the same size in export_file, exporting to file "+file);
	for(size_t i=0; i<v1.size(); i++)
		s << v1[i] << '\t' << v2[i] << '\n';
	s.close();
}

//***************************************************************************
//***************************************************************************