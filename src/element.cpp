#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

#include "element.h"


#define PI      3.1415926535898
#define H_J     6.62607015e-34
#define AMU_KG  1.66053906660e-27
#define KB_J    1.380649e-23


Element element_from_input(const std::string input) {
    std::stringstream ss;
    ss << (input);
    std::string sym;
    double wt, mu, rho, rmin, rmax, padd;
    ss >> sym >> wt >> mu >> rho >> rmin >> rmax >> padd;
    rho /= 2.71828182846;
    return Element(std::move(sym), wt, mu, rho, rmin, rmax, padd);
}

void Element::update_tb(double temperature) {
    tb_ = H_J / sqrt((2*PI*wt_*AMU_KG*KB_J*temperature)) * 1e10;
}

void Element::print() const {
    std::cout << sym_ << '\t' << wt_ << '\t' << tb_ << '\t' << mu_ << '\t';
    std::cout << r_min_ << '\t' << r_max_ << '\t' << p_add_ << std::endl;
}
