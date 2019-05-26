#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include "vec.h"

void vec :: import(double *a)
{
	x[0] = a[0];
	x[1] = a[1];
	x[2] = a[2];
}

void vec :: clean()
{
	x[0] = x[1] = x[2] = 0;
}

vec  vec :: operator+(const vec& B)
{
	vec res;
	res.x[0] = x[0] + B.x[0];
	res.x[1] = x[1] + B.x[1];
	res.x[2] = x[2] + B.x[2];
	return res;
}

vec  vec :: operator-(const vec& B)
{
	vec res;
	res.x[0] = x[0] - B.x[0];
	res.x[1] = x[1] - B.x[1];
	res.x[2] = x[2] - B.x[2];
	return res;
}

vec  vec :: operator*(const double& B)
{
	vec res;
	res.x[0] = x[0]*B;
	res.x[1] = x[1]*B;
	res.x[2] = x[2]*B;
	return res;
}

vec  vec :: operator/(const double& B)
{
	vec res;
	res.x[0] = x[0]/B;
	res.x[1] = x[1]/B;
	res.x[2] = x[2]/B;
	return res;
}

vec & vec :: operator=(const vec& B)
{
	x[0] = B.x[0];
	x[1] = B.x[1];
	x[2] = B.x[2];
	return *this;
}

vec & vec :: operator=(double* B)
{
	x[0] = B[0];
	x[1] = B[1];
	x[2] = B[2];
	return *this;
}

double vec :: operator*(const vec& B)
{
	double res = 0;
	res = x[0]*B.x[0] + x[1]*B.x[1] + x[2]*B.x[2];
	return res;
}

vec vec :: operator^(const vec& B)
{
	vec res;
	res.x[0] = x[1]*B.x[2] - x[2]*B.x[1];
	res.x[1] = x[2]*B.x[0] - x[0]*B.x[2];
	res.x[2] = x[0]*B.x[1] - x[1]*B.x[0];
	return res;
}

std::istream& operator>>(std::istream& in, vec &B)
{
    in>>B.x[0]>>B.x[1]>>B.x[2];
    return in;
}

std::ostream& operator<<(std::ostream& out, vec B)
{
    out<<B.x[0]<<'\t'<<B.x[1]<<'\t'<<B.x[2];
    return out;
}

std::ofstream& operator<<(std::ofstream& out, vec B)
{
    out<<B.x[0]<<'\t'<<B.x[1]<<'\t'<<B.x[2];
    return out;
}

std::stringstream& operator>>(std::stringstream& in, vec &B)
{
    in>>B.x[0]>>B.x[1]>>B.x[2];
    return in;
}

double vec :: norm()
{
	double res=0;
	res = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	return sqrt(res);
}

vec vec :: rand_norm()
{
	x[0] = (double)rand()/RAND_MAX;
	x[1] = (double)rand()/RAND_MAX;
	x[2] = (double)rand()/RAND_MAX;
	while (this->norm() > 1)
	{
		x[0] = (double)rand()/RAND_MAX;
		x[1] = (double)rand()/RAND_MAX;
		x[2] = (double)rand()/RAND_MAX;
	}
	return (*this);
}

//==========debug================

void vec :: print()
{
	std::cout<<x[0]<<'\t'<<x[1]<<'\t'<<x[2]<<std::endl;
}
