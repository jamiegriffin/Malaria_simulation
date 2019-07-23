#include <map>

template<class F>
double find_root(F &f, const double x_guess, const double tol=1E-6, const double r0=2.0, const double r1=1.0);