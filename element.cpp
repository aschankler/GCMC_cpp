#include <iostream>
#include <cstring>
#include <sstream>
#include "element.h"

using namespace std;

void element :: get_param(stringstream& ss)
{

	// get number of element
	/*
	getline(in,tmp);
	while(tmp.find(label_num) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'=');
	ss >> num_el;
	ss.str("");
	ss.clear();
	in.clear();
	in.seekg(ios::beg);
	*/
	// get info of each element
	ss>>sym>>wt>>r_min>>r_max>>p_add;
}

void element :: print()
{
	cout<<sym<<'\t'<<wt<<'\t'<<r_min<<'\t'<<r_max<<'\t'<<p_add<<endl;
}
