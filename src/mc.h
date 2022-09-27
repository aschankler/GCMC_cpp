#ifndef __MC__
#define __MC__

#include <istream>
#include <memory>
#include <vector>
#include "cell.h"


class MCMove {
// Invariant: available must be called before get_new_structure,
//  which must be called before prefactor
  public:
    MCMove(double w, std::string n) : weight_(w), name_(n) {}
    double weight_;
    std::string name_;
    virtual bool available(const Cell&) = 0;
    virtual Cell get_new_structure(const Cell&) = 0;
    virtual double prefactor(const Cell &old, const Cell &trial) = 0;
};


class AddDropMove : public MCMove {
  public:
    AddDropMove(double w = 1.0, std::string n = "AddDrop")
        : MCMove(w, n), prob_add_(0.5), prob_drop_(0.5) {}
    AddDropMove(double wa, double wb, std::string n = "AddDrop")
        : MCMove(wa + wb, n),
          prob_add_(wa / (wa + wb)),
          prob_drop_(wb / (wa + wb)) {}
    bool available(const Cell&) override;
    Cell get_new_structure(const Cell&) override;
    double prefactor(const Cell &old, const Cell &trial) override;
  private:
    double prob_add_;
    double prob_drop_;
    bool add_avail_ = false;
    bool drop_avail_ = false;
    // Cache for add move
    int add_type_;
    vec add_pos_;
    // Cache for drop move
    int drop_type_;
};


class SwapMove : public MCMove {
  public:
    SwapMove(double w = 1., std::string n = "Swap") : MCMove(w, n) {}
    bool available(const Cell&) override;
    Cell get_new_structure(const Cell&) override;
    double prefactor(const Cell &old, const Cell &trial) override;
};


class DoubleMove : public MCMove {
  public:
    DoubleMove(double w = 1., std::string n = "Double")
        : MCMove(w, n), prob_add_(0.5), prob_drop_(0.5) {}
    DoubleMove(double wa, double wb, std::string n = "Double")
        : MCMove(wa + wb, n),
          prob_add_(wa / (wa + wb)),
          prob_drop_(wb / (wa + wb)) {}
    bool available(const Cell&) override;
    Cell get_new_structure(const Cell&) override;
    double prefactor(const Cell &old, const Cell &trial) override;
  private:
    double prob_add_, prob_drop_;
    bool add_avail_ = false;
    bool drop_avail_ = false;
    // Cache add
    int add_type_[2];
    vec add_pos_[2];
    // Cache drop
    int drop_type_[2];
};


typedef std::vector<std::shared_ptr<MCMove>> MCMoveList;

class MCMC {
  public:
    MCMC(MCMoveList moves, double temp, int max_it, bool test = false)
        : moves_(moves), temperature(temp), max_iter(max_it), if_test(test) {}
    // maximum number of iteration
    int max_iter;
    // simulation temperature
    double temperature;
    // if runing test
    bool if_test;

    // global optimized formation energy
    double opt_e = 0.;
    // global optimized cell
    Cell opt_c;

    std::shared_ptr<MCMove> choose_next_move(const Cell&);
    std::shared_ptr<MCMove> get_last_move() const;

    Cell create_new_structure(const Cell&);
    void save_opt_structure(const Cell&);
    bool check_if_accept(Cell& c_old, Cell& c_new);

    void log_header(std::ostream&, const Cell&) const;
    void log_step(std::ostream&, const Cell&, int) const;
    void print() const;
  private:
    // All available moves
    MCMoveList moves_;
    // which action is chosen
    int act_type_ = -1;
    // formation energy of old and new structure
    double e_old_ = 0.;
    double e_trial_ = 0.;
    // status of if accepted or not
    bool accept_ = false;
};

MCMC mcmc_from_in(std::istream&);

#endif  // __MC__
