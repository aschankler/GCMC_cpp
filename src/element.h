#ifndef ELEMENT
#define ELEMENT

#include <string>
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
	// reference density
	double rho;
	// min max threshold
	double r_min, r_max;
	// probabilit of choosing to add
	double p_add;
friend atom;
friend cell;
friend mc;
friend qe_cmd;
public:
	void get_param(std::stringstream& ss);
	void update_tb(double T);

	void print();
};

#endif
