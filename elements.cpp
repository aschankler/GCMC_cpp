#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include "elements.h"

using namespace std;

void elements :: get_param(ifstream& in)
{
	string label_num = "num_el";
	string label_el  = "begin_elements";
	string tmp;
	stringstream ss;
	double sum_p;

	// get number of elements
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

	// allocate elements
	sym   = new string[num_el];
	wt    = new double[num_el];
	l_tb  = new double[num_el];
	r_min = new double[num_el];
	r_max = new double[num_el];
	p_add = new double[num_el];

	// get info of each element
	while(tmp.find(label_el) == string::npos)
		getline(in, tmp);
	for(int t1=0; t1<num_el; t1++)
	{
		getline(in,tmp);
		ss << (tmp);
		ss>>sym[t1]>>wt[t1]>>r_min[t1]>>r_max[t1]>>p_add[t1];
		ss.str("");
		ss.clear();
	}

	// normalize p_add
	sum_p = 0;
	for(int t1=0; t1<num_el; t1++)
		sum_p += p_add[t1];
	for(int t1=0; t1<num_el; t1++)
		p_add[t1] /= sum_p;

	// rewind file
	in.clear();
	in.seekg(ios::beg);
}

void elements :: print()
{
	for(int t1=0; t1<num_el; t1++)
		cout<<sym[t1]<<'\t'<<wt[t1]<<'\t'<<r_min[t1]<<'\t'<<r_max[t1]<<'\t'<<p_add[t1]<<endl;
}
