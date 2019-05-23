#include <iostream>
#include <fstream>
#include <cstring>
#include "cell.h"
#include "vec.h"
#include "atom.h"

using namespace std;

void cell :: read_from_in(ifstream& in)
{
	string label_num_ele = "num_ele";
	string label_num_atm = "num_atm";
	string label_ele = "begin_elements";
	string label_lat = "begin_lattice";
	string label_atm = "begin_atom_positions";
	string tmp;
	stringstream ss;

	// get number of elements
	getline(in,tmp);
	while(tmp.find(label_num_ele) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'='); ss >> num_ele;
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);

	// get number of atoms
	getline(in,tmp);
	while(tmp.find(label_num_atm) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'='); ss >> num_atm;
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);

	// generate element list
	ele_list.resize(num_ele);
	getline(in,tmp);
	while(tmp.find(label_ele) == string::npos)
		getline(in,tmp);
	for(size_t t1=0; t1<num_ele; t1++)
	{
		getline(in,tmp);
		ss << (tmp);
		ele_list[t1].get_param(ss);
		ss.str(""); ss.clear();
	}
	in.clear(); in.seekg(ios::beg);

	// generate atom list
	atm_list.resize(num_atm);
	getline(in,tmp);
	while(tmp.find(label_atm) == string::npos)
		getline(in,tmp);
	for(size_t t1=0; t1<num_atm; t1++)
	{
		getline(in,tmp);
		ss << (tmp);
		atm_list[t1].line_from_in(ss,ele_list);
		ss.str(""); ss.clear();
	}
	in.clear(); in.seekg(ios::beg);

	// get lattice parameter
	while(tmp.find(label_lat) == string::npos)
		getline(in,tmp);
	in>>latt[0]>>latt[1]>>latt[2];
	in.clear(); in.seekg(ios::beg);
}

void cell :: print()
{
	cout<<"Number of elements: "<<num_ele<<endl;
	cout<<"Number of atoms: "<<num_atm<<endl;
	cout<<"List of elements:"<<endl;
	for(size_t t1=0; t1<num_ele; t1++)
		ele_list[t1].print();
	cout<<"List of atoms:"<<endl;
	for(size_t t1=0; t1<num_atm; t1++)
		atm_list[t1].print();
}
