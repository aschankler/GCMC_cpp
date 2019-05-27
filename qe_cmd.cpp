#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <sstream>
#include "qe_cmd.h"

using namespace std;

void qe_cmd :: read_from_in(ifstream& in)
{
	string label_num_core = "num_core";
	string label_npool = "npool";
	string label_ndiag = "ndiag";
	string label_qe_exe = "qe_launcher";
	string label_mpi_launcher = "mpi_launcher";
	string tmp;
	stringstream ss;

	// find number of core, npool and ndiag
	getline(in,tmp);
	while(tmp.find(label_num_core) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'='); ss >> num_core;
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);
	//
	while(tmp.find(label_npool) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'='); ss >> npool;
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);
	//
	while(tmp.find(label_ndiag) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'='); ss >> ndiag;
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);
	// find launchers
	while(tmp.find(label_qe_exe) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'"'); getline(ss,qe_exe,'"');
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);
	//
	while(tmp.find(label_mpi_launcher) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'"'); getline(ss,mpi_launcher,'"');
	ss.str(""); ss.clear(); in.clear(); in.seekg(ios::beg);
}

void qe_cmd :: write_qe_in(ifstream& in, ofstream& out, cell& c_new)
{
	string label_qe_input = "begin_qe_input";
	string label_qe_input2 = "end_qe_input";
	string label_nat = "nat";
	string cell_param = "CELL_PARAMETERS (angstrom)";
	string position = "ATOMIC_POSITIONS (angstrom)";
	string tmp;

	getline(in,tmp);
	while(tmp.find(label_qe_input) == string::npos)
		getline(in,tmp);
	getline(in,tmp);
	// copy parameter to qe.in
	while(tmp.find(label_qe_input2) == string::npos)
	{
		if(tmp.find(label_nat) == string::npos)
			out<<tmp<<endl;
		else
			out<<tmp<<" "<<c_new.num_atm<<endl;
		getline(in,tmp);
	}
	in.clear(); in.seekg(ios::beg);
	// copy cell parameter and position to qe.in
	out<<endl<<cell_param<<endl;
	out<<c_new.latt[0]<<endl;
	out<<c_new.latt[1]<<endl;
	out<<c_new.latt[2]<<endl;

	out<<endl<<position<<endl;
	for(size_t t1=0; t1<c_new.num_atm; t1++)
	{
		out<<setw(2)<<c_new.atm_list[t1].ele->sym<<"    ";
		out<<c_new.atm_list[t1].pos;
		if(c_new.atm_list[t1].if_move == 0)
			out<<" 0 0 0"<<endl;
		else
			out<<endl;
	}
}

void qe_cmd :: call(int if_test)
{
	stringstream ss;
	string tmp;
	if (if_test)
		ss<<mpi_launcher<<" -n "<<num_core<<" echo '    Testing'";
	else
		ss<<mpi_launcher<<" -n "<<num_core<<" "<<qe_exe<<" -npool "<<npool<<" -ndiag "<<ndiag<<" -input qe.in > qe.out";
	tmp = ss.str();
	cout<<endl<<"Launching QE using command:"<<endl;
	cout<<"    "<<tmp<<endl;
	system(tmp.c_str());
	cout<<"QE finished"<<endl;
}
