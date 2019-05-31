#ifndef MC
#define MC

#include <fstream>
#include <vector>
#include "class.h"
#include "cell.h"

class mc
{
public:
	// maximum number of iteration
	int max_iter;
	// simulation temperature
	double T;
	// if runing test
	int if_test;
	// action probability (0, add; 1, remove; 2, swap)
	double act_p[3];
	// which action is chosen
	int act_type;
	// number of elements
	int num_ele;
	// number of atoms changed for each element
	vector<int> num_atm_each_change;
	// global optimized formation energy
	double opt_e;
	// global optimized cell
	cell opt_c;
	// formation energy of old and new structure
	double e1,e2;
	// status of if accepted or not
	int accept;
	// parameters
	double ry_ev = 13.605693009;
	double kb = 8.6173303e-5;  // ev/k

	void read_from_in(std::ifstream& in);

	int factor(int n);
	void create_new_structure(cell c_old, cell& c_new);
	void save_opt_structure(cell& c_new);
	int check_if_accept(cell& c_old, cell& c_new);

	void print();
};

#endif
