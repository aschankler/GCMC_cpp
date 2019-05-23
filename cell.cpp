#include <fstream>
#include <cstring>
#include "cell.h"
#include "vec.h"
#include "atom.h"

using namespace std;

void cell :: read_from_in(ifstream& in,elements& el)
{
	string label_num_atm = "num_atm";
	string label_latt = "begin_lattice";
	string label_atom = "begin_atom_position";
	string tmp;
	stringstream ss;
	atom *pp;

	// get number of atoms
	getline(in,tmp);
	while(tmp.find(label_num_atm) == string::npos)
		getline(in,tmp);
	ss << (tmp);
	getline(ss,tmp,'=');
	ss>>num_atm;
	ss.str("");
	ss.clear();
	in.clear();
	in.seekg(ios::beg);

	// get lattice parameter
	while(tmp.find(label_latt) == string::npos)
		getline(in,tmp);
	in>>latt[0]>>latt[1]>>latt[2];
	in.clear();
	in.seekg(ios::beg);

	// get atoms
	while(tmp.find(label_atom) == string::npos)
		getline(in,tmp);
	pp = at_head;
	for(int t1=0; t1<num_atm; t1++)
	{
		pp->next = new atom;
		pp = pp->next;
		getline(in,tmp);
		pp->line_from_in(tmp,el);
	}

	in.clear();
	in.seekg(ios::beg);
}

void cell :: print()
{
	atom *pp = at_head;
	for(int t1=0; t1<num_atm; t1++)
	{
		pp = pp->next;
		pp->print();
	}
}
