#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include "mc.h"

using namespace std;

void mc :: read_from_in(ifstream& in)
{
	string label_max_iter = "max_iter";
	string label_T = "temperature";
	string label_if_test = "if_test";
	string label_act_p = "begin_action_probability";
	string tmp;
	stringstream ss;

	// get max number of iteration
	getline(in, tmp);
	while(tmp.find(label_max_iter) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'='); ss >> max_iter;
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);

	// get running temperature
	getline(in, tmp);
	while(tmp.find(label_T) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'='); ss >> T;
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);

	// get if running test mode (No QE call)
	getline(in, tmp);
	while(tmp.find(label_if_test) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'='); ss >> if_test;
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);

	// get action probability
	getline(in, tmp);
	while(tmp.find(label_act_p) == string::npos)
		getline(in,tmp);
	in>>act_p[0]>>act_p[1]>>act_p[2];
	in.clear(); in.seekg(ios::beg);
}

void mc :: print()
{
	cout<<"Max iteration: "<<max_iter<<endl;
	cout<<"Simulation temperature: "<<T<<endl;
	cout<<"If running test: "<<if_test<<endl;
	cout<<"Action probability: "<<endl;
	cout<<'\t'<<act_p[0]<<'\t'<<act_p[1]<<'\t'<<act_p[2]<<endl;
}
