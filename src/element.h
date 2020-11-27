#ifndef __ELEMENT__
#define __ELEMENT__

#include <string>
#include <sstream>

class element
{
public:
	// element symbol
	std::string sym;
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

	void get_param(std::string tmp);
	void update_tb(double T);

	void print();
};

#endif
