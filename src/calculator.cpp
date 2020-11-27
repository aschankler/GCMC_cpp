#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "calculator.h"
#include "auxiliary.h"

using namespace std;

void calculator :: read_from_in(ifstream& in)
{
	read(in,"num_core",'=',num_core);
	read(in,"npool",'=',npool);
	read(in,"ndiag",'=',ndiag);
	read(in,"qe_launcher",'=',qe_exe);
	read(in,"mpi_launcher",'=',mpi_launcher);
}

void calculator :: write_qe_in(ifstream& in, ofstream& out, cell& c_new)
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
	for(int t1=0; t1<c_new.num_atm; t1++)
	{
		out<<setw(2)<<c_new.atm_list[t1].ele->sym<<"    ";
		out<<c_new.atm_list[t1].pos;
		if(c_new.atm_list[t1].if_move == 0)
			out<<" 0 0 0"<<endl;
		else
			out<<endl;
	}
}

void calculator :: call(int if_test)
{
	stringstream ss;
	string tmp;
	if (if_test)
		ss<<"mpirun"<<" -n "<<num_core<<" echo '    Testing'";
	else
		ss<<"mkdir -p QE_run; mv qe.in QE_run; cp -r *.upf QE_run 2>/dev/null; cd QE_run; "<<mpi_launcher<<" -n "<<num_core<<" "<<qe_exe<<" -npool "<<npool<<" -ndiag "<<ndiag<<" -input qe.in > qe.out; "<<"cp qe.out ../; cd ../";
	tmp = ss.str();
	cout<<endl<<"Launching QE using command:"<<endl;
	cout<<"    "<<tmp<<endl;
	system(tmp.c_str());
	cout<<"QE finished"<<endl;
}
