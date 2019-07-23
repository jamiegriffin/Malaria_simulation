//************************************************************************
//************************************************************************

template<class T>
void copy(const vector<T> &v, vector<T> &w){
	const int m=v.size();
	const int n=w.size();
	if(m != n)
		error_crit("vectors must be the same size in copy()");
	for(int i=0; i<n; i++)
		w[i]=v[i];
}

template <class C>
typename C::value_type min_el(const C &v){
	return *(min_element(v.begin(), v.end()));
}
template <class C>
typename C::value_type max_el(const C &v){
	return *(max_element(v.begin(), v.end()));
}

//************************************************************************
// remove (at most) one element with value x from v 
template<class T>
bool remove_once(vector<T> &v, T &x){
	typename vector<T>::iterator it=find(v.begin(), v.end(), x);
	if(it != v.end()){
		if(it < v.end()-1){
			*it = *(v.end()-1);
		}
		v.pop_back();
		return true;
	}
	return false;
}
#include "import_file.h"

//***************************************************************************
//***************************************************************************
