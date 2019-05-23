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
	size_t num_ele;
	size_t num_atm;
	double energy;
	vec latt[3];
	vector<element> ele_list;
	vector<atom> atm_list;

	void read_from_in(std::ifstream& in);
	void read_from_qe(std::ifstream& in);

	void print();
};


#endif
