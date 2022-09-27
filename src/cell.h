#ifndef __CELL__
#define __CELL__

#include <istream>
#include <vector>

#include "element.h"
#include "atom.h"
#include "vec.h"


class CellControlParams {
  public:
    CellControlParams(double hmin, double hmax)
        : h_min_(hmin), h_max_(hmax) {}
    // threshold of height
    double h_min_, h_max_;
    // threshold for addition in the xy plane
    double a_min_ = 0.0;
    double a_max_ = 1.0;
    double b_min_ = 0.0;
    double b_max_ = 1.0;
    // if vc-relax
    bool if_vc_relax_ = false;
    // if change-v
    bool if_change_v_ = false;
};

CellControlParams cell_control_from_in(std::istream&);


class Cell {
  public:
    Cell(CellControlParams c = CellControlParams(0, 1)) : control(c) {}
    // Calculation parameters
    CellControlParams control;
    // number of elements
    int num_ele = 0;
    // number of atoms
    int num_atm = 0;
    // lattice parameter
    vec latt[3];
    vec latt_inv[3];
    // list of element
    std::vector<Element> ele_list;
    // list of atoms
    std::vector<Atom> atm_list;

    // energy of cell
    double energy = 0;
    // number of movable atoms
    int num_atm_move, num_atm_remove;
    // number of atoms belonging each element
    std::vector<int> num_ele_each;
    // number of movable atoms belonging each element
    std::vector<int> num_ele_each_move, num_ele_each_remove;


    // io related function
    void read_output(int calculator_type);
    friend void read_from_qe(Cell&);
    friend void read_from_vasp(Cell&);

    void write_axsf(std::ofstream& out) const;
    void write_axsf(std::ofstream& out,int iter) const;
    void write_xsf(std::ofstream& out) const;
    void write_xsf(std::ofstream& out,int iter) const;

    // self-opearted functions
    void count_move_atoms();
    double get_volume() const { return vol_; }
    void update_volume();
    void update_lat_inv();
    void zero_force();
    void update_temperature(double);
    const vec to_crystal(const vec& pos) const;
    const vec from_crystal(const vec& pos) const;

    // adjust atom list
    void ad_atom(vec pos, int ele_type);
    void rm_atom(int ind_atm);
    void sp_atom(int s1, int s2);

    // calculation with input
    void update_tb(double T);
    double min_distance(vec) const;
    double min_distance(vec, int&) const;

    void print() const;

  private:
    // volume of cell
    double vol_;
};

Cell cell_from_in(std::istream&);

#endif  // __CELL__
