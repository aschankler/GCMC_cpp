#ifndef __MC__
#define __MC__

#include <fstream>
#include <vector>
#include "cell.h"

class mc {
public:
	// maximum number of iteration
	int max_iter;
	// simulation temperature
	double temperature;
	// if runing test
	int if_test;
	// number of action types
	int const static num_act = 3;
	// action probability (0, add; 1, remove; 2, swap)
	double act_p[num_act];

	int act_type;

	// global optimized formation energy
	double opt_e;
	// global optimized cell
	Cell opt_c;

	void read_from_in(std::ifstream& in);

	int factor(int n);
	void create_new_structure(const Cell c_old, Cell& c_new);
	void save_opt_structure(const Cell c_new);
	int check_if_accept(Cell& c_old, Cell& c_new);

    void log_header(std::ostream&, const Cell&) const;
    void log_step(std::ostream&, const Cell&, int) const;
	void print();
  private:
    // which action is chosen
    int act_type_;
    // formation energy of old and new structure
    double e_old_, e_trial_;
    // number of atoms changed for each element
    std::vector<int> num_atm_each_change_;
    // status of if accepted or not
    bool accept_;
};

#endif  // __MC__
