#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include "cell.h"
#include "vec.h"
#include "atom.h"

using namespace std;

void cell :: read_from_in(ifstream& in)
{
	string label_num_ele = "num_ele";
	string label_num_atm = "num_atm";
	string label_h_min = "h_min";
	string label_h_max = "h_max";
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

	// get threshold of height
	getline(in,tmp);
	while(tmp.find(label_h_min) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'='); ss >> h_min;
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);

	getline(in,tmp);
	while(tmp.find(label_h_max) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'='); ss >> h_max;
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

void cell :: read_from_qe(ifstream& in)
{
	string label_final_energy = "Final energy";
	string label_final_position = "Begin final coordinates";
	string tmp;
	stringstream ss;
	getline(in,tmp);
	// find final status
	while(tmp.find(label_final_energy) == string::npos && !in.eof())
		getline(in,tmp);
	if (! in.eof())
	{
		// save final energy
		ss << (tmp);
		getline(ss,tmp,'=');
		ss >> energy;
		ss.str(""); ss.clear();
		// save final postion
		while(tmp.find(label_final_position) == string::npos)
			getline(in,tmp);
		getline(in,tmp);
		getline(in,tmp);
		for (size_t t1=0; t1<num_atm; t1++)
		{
			in>>tmp>>atm_list[t1].pos;
			getline(in,tmp);
		}
	}
	else
	{
		cout<<"Warning: QE calculation does not converge, set energy to 0 Ry"<<endl;
		energy = 0;
	}
}

void cell :: count_move_atoms()
{
	num_ele_each.resize(num_ele);
	num_ele_each_move.resize(num_ele);
	num_atm_move = 0;
	for(size_t t1=0; t1<num_ele; t1++)
	{
		num_ele_each[t1] = 0;
		num_ele_each_move[t1] = 0;
	}
	for(size_t t1=0; t1<num_atm; t1++)
	{
		num_atm_move += atm_list[t1].if_move;
		num_ele_each[atm_list[t1].type]++;
		num_ele_each_move[atm_list[t1].type] += atm_list[t1].if_move;
	}
}

double cell :: get_volume()
{
	double h[3]={0,0,h_max-h_min};
	vec hh;
	hh=&h[0];
	vol = (latt[0]^latt[1])*hh;
	return vol;
}

void cell :: ad_atom(vec pos, int ele_type)
{
	atom tmp;
	tmp.type = ele_type;
	tmp.ele = &ele_list[ele_type];
	tmp.pos = pos;
	tmp.force = pos*0;
	tmp.if_move = 1;
	atm_list.push_back(tmp);
	num_atm++;
	count_move_atoms();
}

void cell :: rm_atom(int ind_atm)
{
	atm_list.erase(atm_list.begin() + ind_atm);
	num_atm--;
	count_move_atoms();
}

void cell :: sp_atom(int s1, int s2)
{
	vec tmp;
	tmp = atm_list[s1].pos;
	atm_list[s1].pos = atm_list[s2].pos;
	atm_list[s2].pos = tmp;
}

void cell :: update_tb(double T)
{
	for(size_t t1=0; t1<num_ele; t1++)
		ele_list[t1].update_tb(T);
}

void cell :: min_distance(vec pos, double& rr, int& ind)
{
	double r_tmp;
	rr = 1e10;
	for(int t1=-1; t1<2; t1++)
	for(int t2=-1; t2<2; t2++)
	for(int t3=-1; t3<2; t3++)
		for(size_t t4=0; t4<num_atm; t4++)
		{
			r_tmp=(pos + latt[0]*t1 + latt[1]*t2 + latt[2]*t3 - atm_list[t4].pos).norm();
			if (rr > r_tmp)
			{
				rr = r_tmp;
				ind = t4;
			}
		}
}

void cell :: print()
{
	cout<<"Number of elements: "<<num_ele<<endl;
	cout<<"Number of atoms: "<<num_atm<<endl;
	cout<<"Threshold of height: "<<h_min<<", "<<h_max<<endl;
	cout<<"Energy: "<<energy<<endl;
	cout<<"Volume: "<<vol<<endl;
	cout<<"Total movable atoms: "<<num_atm_move<<endl;
	cout<<"Atoms per elements: "<<endl;
	for(size_t t1=0; t1<num_ele; t1++)
		cout<<ele_list[t1].sym<<'\t'<<num_ele_each[t1]<<endl;
	cout<<"Movable atoms per elements: "<<endl;
	for(size_t t1=0; t1<num_ele; t1++)
		cout<<ele_list[t1].sym<<'\t'<<num_ele_each_move[t1]<<endl;
//	cout<<"List of elements:"<<endl;
//	for(size_t t1=0; t1<num_ele; t1++)
//		ele_list[t1].print();
	cout<<"List of atoms:"<<endl;
	for(size_t t1=0; t1<num_atm; t1++)
		atm_list[t1].print();
}
