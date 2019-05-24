#include <iostream>
#include <cstring>
#include <sstream>
#include "element.h"

using namespace std;

void element :: get_param(stringstream& ss)
{
	ss>>sym>>wt>>mu>>r_min>>r_max>>p_add;
}

void element :: print()
{
	cout<<sym<<'\t'<<wt<<'\t'<<mu<<'\t'<<r_min<<'\t'<<r_max<<'\t'<<p_add<<endl;
}
