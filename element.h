#ifndef ELEMENT
#define ELEMENT

#include <cstring>
#include <sstream>
#include "class.h"

using namespace std;

class element
{
private:
	// element symbol
	string sym;
	// atomic weight
	double wt;
	// thermal debroye wavelength
	double l_tb;
	// min max threshold
	double r_min, r_max;
	// probabilit of choosing to add
	double p_add;
friend atom;
friend cell;
public:
	void get_param(std::stringstream& ss);

	void print();
};

#endif
