#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include "vec.h"
#include "rng.h"

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

vec  vec :: operator+(const vec& B) const
{
	vec res;
	res.x[0] = x[0] + B.x[0];
	res.x[1] = x[1] + B.x[1];
	res.x[2] = x[2] + B.x[2];
	return res;
}

vec  vec :: operator-(const vec& B) const
{
	vec res;
	res.x[0] = x[0] - B.x[0];
	res.x[1] = x[1] - B.x[1];
	res.x[2] = x[2] - B.x[2];
	return res;
}

vec  vec :: operator*(const double& B) const
{
	vec res;
	res.x[0] = x[0]*B;
	res.x[1] = x[1]*B;
	res.x[2] = x[2]*B;
	return res;
}

vec  vec :: operator/(const double& B) const
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

double vec :: operator*(const vec& B) const
{
	double res = 0;
	res = x[0]*B.x[0] + x[1]*B.x[1] + x[2]*B.x[2];
	return res;
}

vec vec :: operator^(const vec& B) const
{
	vec res;
	res.x[0] = x[1]*B.x[2] - x[2]*B.x[1];
	res.x[1] = x[2]*B.x[0] - x[0]*B.x[2];
	res.x[2] = x[0]*B.x[1] - x[1]*B.x[0];
	return res;
}

const double vec::operator[](const std::size_t idx) const {
    return x[idx];
}

double& vec::operator[](const std::size_t idx) {
    return x[idx];
}

std::istream& operator>>(std::istream& in, vec &B)
{
    in>>B.x[0]>>B.x[1]>>B.x[2];
    return in;
}

std::ostream& operator<<(std::ostream& out, const vec &B)
{
    out<<B.x[0]<<'\t'<<B.x[1]<<'\t'<<B.x[2];
    return out;
}

std::ofstream& operator<<(std::ofstream& out, const vec &B)
{
    //out<<B.x[0]<<'\t'<<B.x[1]<<'\t'<<B.x[2];
    out<<std::fixed<<std::setw(22)<<std::setprecision(15)<<B.x[0]<<std::setw(22)<<std::setprecision(15)<<B.x[1]<<std::setw(22)<<std::setprecision(15)<<B.x[2];
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


vec &vec::rand() {
    x[0] = gcmc::rand_uniform();
    x[1] = gcmc::rand_uniform();
    x[2] = gcmc::rand_uniform();
    return *this;
}


vec &vec :: rand_norm()
{
    this->rand();
    while (this->norm() > 1) {
        this->rand();
    }
    return *this;
}

//==========debug================

void vec :: print() const {
	std::cout<<x[0]<<'\t'<<x[1]<<'\t'<<x[2]<<std::endl;
}
