#ifndef CELL
#define CELL

#include <fstream>
#include <vector>
#include "class.h"
#include "vec.h"
#include "element.h"
#include "atom.h"

class cell
{
public:
	int num_ele;
	int num_atm;
	double energy;
	vec latt[3];
	vector<element> ele_list;
	vector<atom> atm_list;
	int num_ele_move;
	vector<int> num_ele_each_move;

	// function relying on external files
	void read_from_in(std::ifstream& in);
	void read_from_qe(std::ifstream& in);

	// self-opearted functions
	void count_move_atoms();

	// adjust atom list
	void ad_atom(vec pos, int ele_type);
	void rm_atom(int ind_atm);
	void sp_atom(int s1, int s2);

	void print();
};


#endif
