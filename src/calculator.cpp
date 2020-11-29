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
	read(in,"exe",'=',exe);
	if (exe.find("pw") != string::npos)
		calculator_type = 1;
	else if(exe.find("vasp") != string::npos)
		calculator_type = 2;
	else
	{
		cout<<"Error: Currently only support QE and VASP!"<<endl;
		exit(EXIT_FAILURE);
	}
	read(in,"mpi_launcher",'=',mpi_launcher);
	read(in,"num_core",'=',num_core);
	switch(calculator_type)
	{
	case 1:
	{
		read(in,"npool",'=',npool);
		read(in,"ndiag",'=',ndiag);
		break;
	}
	case 2:
	{
		break;
	}
	default:
	{
		cout<<"Error: Invalid calculator type!"<<endl;
		exit(EXIT_FAILURE);
	}
	}
}

void calculator :: write_input(ifstream& in, cell& c_new)
{
	switch (calculator_type)
	{
		case 1: //QE
		{
			write_qe_in(in,c_new);
			break;
		}
		case 2: //vasp
		{
			write_vasp_in(in,c_new);
			break;
		}
		default:
		{
			cout<<"Error: Invalid calculator type!"<<endl;
			exit(EXIT_FAILURE);
		}
	}
}

void calculator :: write_qe_in(ifstream& in, cell& c_new)
{
	ofstream out("qe.in");
	string label_qe_input = "begin_qe_input";
	string label_qe_input2 = "end_qe_input";
	string label_nat = "nat";
	string cell_param = "CELL_PARAMETERS (angstrom)";
	string position = "ATOMIC_POSITIONS (angstrom)";
	string tmp;

	in.clear(); in.seekg(ios::beg);
	while(getline(in,tmp))
		if(tmp.find(label_qe_input) != string::npos)
			break;
	if(in.eof())
	{
		cout<<"Error: Cannot find "<<label_qe_input<<'!'<<endl;
		exit(EXIT_FAILURE);
	}
	// copy parameter to qe.in
	while(getline(in,tmp))
	{
		if(tmp.find(label_qe_input2) == string::npos)
		{
			if(tmp.find(label_nat) == string::npos)
				out<<tmp<<endl;
			else
				out<<tmp<<" "<<c_new.num_atm<<endl;
		}
		else
			break;
	}
	if(in.eof())
	{
		cout<<"Error: Cannot find "<<label_qe_input2<<'!'<<endl;
		exit(EXIT_FAILURE);
	}
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
	out.close();
}

void calculator ::write_vasp_in(std::ifstream &in, cell &c_new)
{
	ofstream incar("INCAR");
	ofstream kpoints("KPOINTS");
	ofstream poscar("POSCAR");
	string label_vasp_incar = "begin_vasp_incar";
	string label_vasp_incar2= "end_vasp_incar";
	string label_vasp_kpoints = "begin_vasp_kpoints";
	string label_vasp_kpoints2= "end_vasp_kpoints";
	string tmp;
	// generate incar
	in.clear(); in.seekg(ios::beg);
	while(getline(in,tmp))
		if(tmp.find(label_vasp_incar) != string::npos)
			break;
	if(in.eof())
	{
		cout<<"Error: Cannot find "<<label_vasp_incar<<'!'<<endl;
		exit(EXIT_FAILURE);
	}
	while(getline(in,tmp))
	{
		if(tmp.find(label_vasp_incar2) == string::npos)
			incar<<tmp<<endl;
		else
			break;
	}
	if(in.eof())
	{
		cout<<"Error: Cannot find "<<label_vasp_incar2<<'!'<<endl;
		exit(EXIT_FAILURE);
	}
	incar.close();
	// generate kpoint
	in.clear(); in.seekg(ios::beg);
	while(getline(in,tmp))
		if(tmp.find(label_vasp_kpoints) != string::npos)
			break;
	if(in.eof())
	{
		cout<<"Error: Cannot find "<<label_vasp_kpoints<<'!'<<endl;
		exit(EXIT_FAILURE);
	}
	while(getline(in,tmp))
	{
		if(tmp.find(label_vasp_kpoints2) == string::npos)
			kpoints<<tmp<<endl;
		else
			break;
	}
	if(in.eof())
	{
		cout<<"Error: Cannot find "<<label_vasp_kpoints2<<'!'<<endl;
		exit(EXIT_FAILURE);
	}
	kpoints.close();

	// generate poscar
	poscar<<"gcmc_auto_generated"<<endl;
	poscar<<"   1.000000000000"<<endl;
	poscar<<c_new.latt[0]<<endl;
	poscar<<c_new.latt[1]<<endl;
	poscar<<c_new.latt[2]<<endl;
	for(auto m1:c_new.ele_list)
		poscar<<setw(5)<<m1.sym;
	poscar<<endl;
	for(auto m1:c_new.num_ele_each)
		poscar<<setw(5)<<m1;
	poscar<<endl;
	poscar<<"Selective dynamics"<<endl;
	poscar<<"Direct"<<endl;
	// atomic positions
	vec frac_coord;
	for(auto m1:c_new.ele_list)
	{
		for(auto m2:c_new.atm_list)
		{
			if(m1.sym == m2.ele->sym)
			{
				frac_coord.x[0] = m2.pos*c_new.latt_inv[0];
				frac_coord.x[1] = m2.pos*c_new.latt_inv[1];
				frac_coord.x[2] = m2.pos*c_new.latt_inv[2];
				poscar<<frac_coord;
				if(m2.if_move)
					poscar<<"  T  T  T"<<endl;
				else
					poscar<<"  F  F  F"<<endl;
			}
		}
	}
	poscar.close();
}

void calculator :: call(int if_test)
{
	stringstream ss;
	string tmp;
	if (if_test)
		ss<<"mpirun"<<" -n "<<num_core<<" echo '    Testing'";
	else
	{
		switch(calculator_type)
		{
		case 1: //QE
		{
			ss<<mpi_launcher<<" -n "<<num_core<<" "<<exe<<" -npool "<<npool<<" -ndiag "<<ndiag<<" -input qe.in > qe.out";
			tmp = ss.str();
			cout<<endl<<"Launching QE using command:"<<endl;
			cout<<"    "<<tmp<<endl;
			system(tmp.c_str());
			cout<<"QE finished"<<endl;
			break;
		}
		case 2: //vasp
		{
			ss<<mpi_launcher<<" -n "<<num_core<<" "<<exe<<" > vasp.out";
			tmp = ss.str();
			cout<<endl<<"Launching VASP using command:"<<endl;
			cout<<"    "<<tmp<<endl;
			system(tmp.c_str());
			cout<<"VASP finished"<<endl;
			break;
		}
		default:
		{
			cout<<"Error: Invalid calculator type!"<<endl;
			exit(EXIT_FAILURE);
		}
		}
	}
}
