#ifndef MY_VEC
#define MY_VEC

#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include "class.h"

class vec
{
private:
	double x[3];
friend atom;
friend cell;
friend mc;
public:
	void import(double *);
	void clean();					// reset value to zero
	vec operator+(const vec&);
	vec operator-(const vec&);
	vec operator*(const double&);
	vec operator/(const double&);
	double operator*(const vec&);
	vec operator^(const vec&);		// for cross product
	vec & operator=(const vec&);
	vec & operator=(double*);		// can replace import
	friend std::istream& operator>>(std::istream&,vec&);
	friend std::ostream& operator<<(std::ostream&,vec);
	friend std::ofstream& operator<<(std::ofstream&,vec);
	friend std::stringstream& operator>>(std::stringstream&,vec&);
	double norm();
	vec rand_norm();

	//debug
	void print();
};

#endif
