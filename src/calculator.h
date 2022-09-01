#ifndef __CALCULATOR__
#define __CALCULATOR__

#include <fstream>
#include <string>
#include "cell.h"

class Calculator {
  public:
    // Needed for encapsulation; should eventually be factored out entirely
    int get_type() const { return calculator_type_; }
    std::string cli_string() const;

    void read_from_in(std::ifstream& in);
    void write_input(std::ifstream& in, const Cell& c_new) const;
    void write_qe_in(std::ifstream& in, const Cell& c_new) const;
    void write_vasp_in(std::ifstream& in, const Cell& c_new) const;
    void call(int if_test) const;

  protected:
    int num_core_;
    int calculator_type_; // 1: qe, 2: vasp
    std::string exe_name_;
    std::string mpi_launcher_;
    std::string para_params_;
};

#endif  // __CALCULATOR__
