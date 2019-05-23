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
private:
	size_t num_ele;
	size_t num_atm;
	vec latt[3];
	vector<element> ele_list;
	vector<atom> atm_list;
public:
	void read_from_in(std::ifstream& in);

	void print();
};


#endif
