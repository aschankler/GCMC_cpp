#include <iostream>
#include <cstring>
#include <sstream>
#include <cmath>
#include "element.h"

using namespace std;

void element :: get_param(stringstream& ss)
{
	ss>>sym>>wt>>mu>>rho>>r_min>>r_max>>p_add;
}

void element :: update_tb(double T)
{
	double h = 6.626070040e-34;
	double kb = 1.38064852e-23;
	double amu_kg = 1.660539040e-27;
	double pi = 3.1415926535898;
	tb=h / sqrt((2*pi*wt*amu_kg*kb*T)) * 1e10;
}

void element :: print()
{
	cout<<sym<<'\t'<<wt<<'\t'<<tb<<'\t'<<mu<<'\t'<<r_min<<'\t'<<r_max<<'\t'<<p_add<<endl;
}
