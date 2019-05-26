#ifndef MC
#define MC

#include <fstream>
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
	// number of atoms changed (1, add; -1, remove; 0, swap)
	int num_change;
	// element type of atoms changed
	int el_change;

	void read_from_in(std::ifstream& in);

	void create_new_structure(cell c_old, cell& c_new);

	void print();
};

#endif
