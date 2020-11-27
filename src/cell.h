#ifndef __CELL__
#define __CELL__

#include <fstream>
#include <vector>
#include "element.h"
#include "atom.h"
#include "vec.h"

class cell
{
public:
	// number of elements
	int num_ele;
	// number of atoms
	int num_atm;
	// lattice parameter
	vec latt[3];
	// list of element
	std::vector<element> ele_list;
	// list of atoms
	std::vector<atom> atm_list;
	// threshold of height
	double h_min, h_max;
	// energy of cell
	double energy;
	// number of movable atoms
	int num_atm_move, num_atm_remove;
	// number of atoms belonging each element
	std::vector<int> num_ele_each;
	// number of movable atoms belonging each element
	std::vector<int> num_ele_each_move, num_ele_each_remove;
	// volume of cell
	double vol;
	// if vc-relax
	int if_vc_relax;
	// if change-v
	int if_change_v;
	// parameters
	double ry_ev = 13.605693009;

	// io related function
	void read_from_in(std::ifstream& in);
	void read_from_qe(std::ifstream& in);
	void write_axsf(std::ofstream& out);
	void write_axsf(std::ofstream& out,int iter);
	void write_xsf(std::ofstream& out);
	void write_xsf(std::ofstream& out,int iter);

	// self-opearted functions
	void count_move_atoms();
	double get_volume();

	// adjust atom list
	void ad_atom(vec pos, int ele_type);
	void rm_atom(int ind_atm);
	void sp_atom(int s1, int s2);

	// calculation with input
	void update_tb(double T);
	void min_distance(vec pos, double& rr, int& ind);

	void print();
};

#endif
