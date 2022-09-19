#ifndef __MY_VEC__
#define __MY_VEC__

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

class vec
{
public:
	double x[3];

	void import(double *);
	void clean();					// reset value to zero
	vec operator+(const vec&) const;
	vec operator-(const vec&) const;
	vec operator*(const double&) const;
	vec operator/(const double&) const;
	double operator*(const vec&) const;
	vec operator^(const vec&) const;  // for cross product
    // Subscripts get at vector components
    double& operator[](const std::size_t);
    const double operator[](const std::size_t) const;
	vec & operator=(const vec&);
	vec & operator=(double*);		// can replace import
	friend std::istream& operator>>(std::istream&, vec&);
	friend std::ostream& operator<<(std::ostream&, const vec&);
	friend std::ofstream& operator<<(std::ofstream&, const vec&);
	friend std::stringstream& operator>>(std::stringstream&, vec&);
	double norm();
	vec &rand();
	vec &rand_norm();

	//debug
	void print() const;
};

#endif
