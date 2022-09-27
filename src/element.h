#ifndef __ELEMENT__
#define __ELEMENT__

#include <string>

class Element {
  public:
    Element(std::string sym, double wt, double mu, double rho, double rmin, double rmax, double p_add)
        : sym_(sym), wt_(wt), mu_(mu), rho_(rho), r_min_(rmin), r_max_(rmax), p_add_(p_add) {}
    // element symbol
    std::string sym_;
    // atomic weight (amu)
    double wt_;
    // chemical potential (eV)
    double mu_;
    // thermal debroye wavelength (Angstrom)
    double tb_;
    // reference density
    double rho_;
    // min max threshold (Angstrom)
    double r_min_, r_max_;
    // probabilit of choosing to add
    double p_add_;

    void update_tb(double);

    void print() const;
};

Element element_from_input(const std::string);

#endif
