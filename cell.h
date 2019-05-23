#ifndef CELL
#define CELL

#include <fstream>
#include "class.h"
#include "vec.h"
#include "atom.h"

class cell
{
private:
	int num_atm;
	vec latt[3];
	atom *at_head;
public:
	void read_from_in(std::ifstream& in, elements& el);

	void print();
};


#endif
