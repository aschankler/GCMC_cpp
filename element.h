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
	// chemical potential
	double mu;
	// thermal debroye wavelength
	double tb;
	// min max threshold
	double r_min, r_max;
	// probabilit of choosing to add
	double p_add;
friend atom;
friend cell;
friend mc;
public:
	void get_param(std::stringstream& ss);
	void update_tb(double T);

	void print();
};

#endif
