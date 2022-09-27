#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <random>
#include <map>
#include "mc.h"
#include "auxiliary.h"
#include "rng.h"


#define EPS  1e-8
// In eV/K
#define KB_EV 8.6173303e-5

using namespace std;

// ----------------------------------------------------------------
// AddDropMove implementation
// ----------------------------------------------------------------

static bool choose_add(const Cell &cell_old, vec &pos_add, int &ele_type) {
    const int num_trial = 1000;

    // Choose species to add
    double add_weight_sum = 0;
    for (int t1 = 0; t1 < cell_old.num_ele; t1++)
        add_weight_sum += cell_old.ele_list[t1].p_add_;

    // Choose random number from [0, add_weight_sum)
    ele_type = cell_old.num_ele;
    std::uniform_real_distribution<double> ele_dist(0., add_weight_sum);
    add_weight_sum = ele_dist(gcmc::rng);
    for (int t1 = 0; t1 < cell_old.num_ele; t1++) {
        add_weight_sum -= cell_old.ele_list[t1].p_add_;
        if(add_weight_sum <= 0) {
            ele_type = t1;
            break;
        }
    }

    Element ele_tmp = cell_old.ele_list[ele_type];

    // Random distribution for new atom addition
    std::uniform_real_distribution<double> da(cell_old.a_min, cell_old.a_max);
    std::uniform_real_distribution<double> db(cell_old.b_min, cell_old.b_max);
    // Breaks for "non standard" lattice orrientations
    double cmin = cell_old.to_crystal({0, 0, cell_old.h_min})[2];
    double cmax = cell_old.to_crystal({0, 0, cell_old.h_max})[2];
    std::uniform_real_distribution<double> dc(cmin, cmax);

    // Choose position to add
    for (int t1 = 0; t1 < num_trial; t1++) {
        // Generate new position
        pos_add = {da(gcmc::rng), db(gcmc::rng), dc(gcmc::rng)};

        // Old atom-centered version
        // pos_add = atm_tmp.pos + pos_add.rand_norm()*(ele_tmp.r_min + (double)rand()/RAND_MAX*(ele_tmp.r_max - ele_tmp.r_min));

        // Filter xy-coordinate range
        if (not (cell_old.a_min < pos_add[0] && pos_add[0] < cell_old.a_max
                 && cell_old.b_min < pos_add[1] && pos_add[1] < cell_old.b_max)
           ) {
            continue;
        }

        // Convert to cartesian coordinates and check Z coordinate
        pos_add = cell_old.from_crystal(pos_add);
        if (pos_add.x[2] < cell_old.h_min || pos_add.x[2] > cell_old.h_max) {
            continue;
        }

        // Filter coordination rule
        if (!cell_old.if_vc_relax) {
            double r_tmp = cell_old.min_distance(pos_add);
            if (r_tmp > ele_tmp.r_min_ && r_tmp < ele_tmp.r_max_) {
                // Accept this position
                return true;
            }
        } else {
            return true;
        }
    }

    // Did not find a suitable posititon
    return false;
}


bool AddDropMove::available(const Cell &cell) {
    // Drop available
    drop_avail_ = (cell.num_atm_remove > 0);

    // Check if add available and cache position
    add_avail_ = choose_add(cell, add_pos_, add_type_);

    return drop_avail_ || add_avail_;
}


static int choose_drop(const Cell &cell_old, int &rm_idx) {
    std::uniform_int_distribution<> dist(0, cell_old.num_atm_remove - 1);
    // Select the ith removable atom
    int atm_id_tmp = dist(gcmc::rng);
    for (int t1 = 0; t1 < cell_old.num_atm; t1++) {
        if (cell_old.atm_list[t1].if_move == 2) {
            // Atom is removable
            if (atm_id_tmp == 0) {
                rm_idx = t1;
                return 1;
            }
            atm_id_tmp--;
        }
    }
    return 0;
}


Cell AddDropMove::get_new_structure(const Cell &cell) {
    Cell cell_new = cell;

    if (add_avail_ and drop_avail_) {
        // Choose to add or drop
        std::uniform_real_distribution<double> dist(0., prob_add_ + prob_drop_);
        if (dist(gcmc::rng) < prob_add_) {
            drop_avail_ = false;  // Disable drop move
        } else {
            add_avail_ = false;
        }
    }

    if (add_avail_) {
        // Use cached addition
        cell_new.ad_atom(add_pos_, add_type_);
        cout << "Atom added, +" << cell_new.num_atm;
        cout << " " << cell_new.ele_list[add_type_].sym_;
        cout << ", position is " << add_pos_ << endl;
        return cell_new;
    }

    if (drop_avail_) {
        int rm_idx;
        choose_drop(cell_new, rm_idx);
        cell_new.rm_atom(rm_idx);
        // Cache move
        drop_type_ = cell.atm_list[rm_idx].type;
        cout << "Atom removed, -" << rm_idx+1;
        cout << " " << cell.atm_list[rm_idx].ele->sym_ << endl;
        return cell_new;
    }

    throw std::runtime_error("Neither add nor drop valid");
}


static double counting_factor(int a, int b) {
    // a! / b!
    if (a < 0 || b < 0) {
        throw std::runtime_error("counting factor can not be negative");
    }

    if (a == b) { return 1; }

    if (a > b) {
        int prod = 1;
        for (int n = a; n > b; n--) {
            prod *= n;
        }
        return prod;
    } else {
        // a < b
        return 1. / counting_factor(b, a);
    }
}

static double phase_space_volume(const Cell &old, const Cell &trial) {
    double phase_vol = 1.;
    for (int i = 0; i < old.num_ele; i++) {
        auto e1 = old.ele_list[i];
        int dn = trial.num_ele_each[i] - old.num_ele_each[i];
        // Factor of (V/Lambda^3)^dN
        phase_vol *= pow(trial.get_volume() / pow(e1.tb_, 3), dn);
        // Counting factor for indistinguishable atoms
        phase_vol *= counting_factor(old.num_ele_each[i], trial.num_ele_each[i]);
    }
    // if allow V to change, multiply by (V_new/V_old)^N_tot(old)
    if(trial.if_change_v)
        phase_vol *= pow((trial.get_volume() / old.get_volume()), old.num_atm);

    return phase_vol;
}


double AddDropMove::prefactor(const Cell &old, const Cell &trial) {
    double prefactor = phase_space_volume(old, trial);

    double add_weight_tot = 0.;
    for (auto elem : trial.ele_list)
        add_weight_tot += elem.p_add_;

    // Prob of choosing reverse move over forward move
    if (add_avail_) {
        // Forward move
        prefactor /= prob_add_ * trial.ele_list[add_type_].p_add_ / add_weight_tot;
        // Reverse move
        prefactor *= prob_drop_ * trial.num_ele_each_remove[add_type_] / trial.num_atm_remove;
    } else if (drop_avail_) {
        // Forward move
        prefactor /= prob_drop_ * old.num_ele_each_remove[drop_type_] / old.num_atm_remove;
        // Reverse move
        prefactor *= prob_add_ * old.ele_list[drop_type_].p_add_ / add_weight_tot;
    } else {
        throw std::runtime_error("Previous move not recorded");
    }

    return prefactor;
}


// ----------------------------------------------------------------
// SwapMove implementation
// ----------------------------------------------------------------

bool SwapMove::available(const Cell &cell) {
    int movable_types = 0;
    for (int i = 0; i < cell.num_ele; i++) {
        if (cell.num_ele_each_move[i] > 0)
            movable_types++;
    }
    return movable_types >= 2;
}


Cell SwapMove::get_new_structure(const Cell &cell) {
    // find the first atom to switch
    std::uniform_int_distribution<> dist(0, cell.num_atm_move-1);
    int idx1 = -1;
    int atm_id_tmp = dist(gcmc::rng);
    for (int t1 = 0; t1 < cell.num_atm; t1++) {
        if (cell.atm_list[t1].if_move >= 1) {
            // Atom is movable
            if (atm_id_tmp == 0) {
                idx1 = t1;
                break;
            }
            atm_id_tmp--;
        }
    }

    if (idx1 < 0)
        throw std::runtime_error("Could not find swap move");

    // find the second atom to switch
    int idx2 = -1;
    int num_move_distinct = cell.num_atm_move - cell.num_ele_each_move[cell.atm_list[idx1].type];
    atm_id_tmp = dist(gcmc::rng, decltype(dist)::param_type(0, num_move_distinct-1));
    for (int t1 = 0; t1 < cell.num_atm; t1++) {
        if (cell.atm_list[t1].if_move >= 1
            && cell.atm_list[t1].type != cell.atm_list[idx1].type
        ) {
            // Atom is movable
            if (atm_id_tmp == 0) {
                idx2 = t1;
                break;
            }
            atm_id_tmp--;
        }
    }

    if (idx2 < 0)
        throw std::runtime_error("Could not find swap move");

    // Create new cell
    cout << "Atoms swapped, " << idx1+1 << " " << cell.atm_list[idx1].ele->sym_;
    cout << " <--> " << idx2+1 << " " << cell.atm_list[idx2].ele->sym_ << endl;
    Cell cell_new = cell;
    cell_new.sp_atom(idx1, idx2);
    return cell_new;
}


double SwapMove::prefactor(const Cell &old, const Cell &trial) {
    // Only volume change affects sampling dist
    double prefactor = 1.;

    if (trial.if_change_v)
        prefactor *= pow((trial.get_volume() / old.get_volume()), old.num_atm);

    return prefactor;
}


// ----------------------------------------------------------------
// DoubleMove implementation
// ----------------------------------------------------------------

bool DoubleMove::available(const Cell &cell) {
    // Drop needs two removable atoms
    drop_avail_ = (cell.num_atm_remove >=2);

    // Check if add available and cache position
    add_avail_ = false;
    if (choose_add(cell, add_pos_[0], add_type_[0])) {
        // First add is available
        Cell cell_tmp = cell;
        cell_tmp.ad_atom(add_pos_[0], add_type_[0]);
        // Choose second add
        add_avail_ = choose_add(cell_tmp, add_pos_[1], add_type_[1]);
    }

    return drop_avail_ || add_avail_;
}


Cell DoubleMove::get_new_structure(const Cell &cell) {
    Cell cell_new = cell;

    if (add_avail_ and drop_avail_) {
        // Choose to add or drop
        std::uniform_real_distribution<double> dist(0., prob_add_ + prob_drop_);
        if (dist(gcmc::rng) < prob_add_) {
            drop_avail_ = false;  // Disable drop move
        } else {
            add_avail_ = false;
        }
    }

    if (add_avail_) {
        // Use cached addition
        cell_new.ad_atom(add_pos_[0], add_type_[0]);
        cout << "Atom added, +" << cell_new.num_atm;
        cout << " " << cell_new.ele_list[add_type_[0]].sym_;
        cout << ", position is " << add_pos_[0] << endl;
        cell_new.ad_atom(add_pos_[1], add_type_[1]);
        cout << "Atom added, +" << cell_new.num_atm;
        cout << " " << cell_new.ele_list[add_type_[1]].sym_;
        cout << ", position is " << add_pos_[1] << endl;
        return cell_new;
    }

    if (drop_avail_) {
        int rm_idx1, rm_idx2;
        choose_drop(cell_new, rm_idx1);
        drop_type_[0] = cell_new.atm_list[rm_idx1].type;
        cout << "Atom removed, -" << rm_idx1+1;
        cout << " " << cell_new.atm_list[rm_idx1].ele->sym_ << endl;
        cell_new.rm_atom(rm_idx1);
        choose_drop(cell_new, rm_idx2);
        drop_type_[1] = cell_new.atm_list[rm_idx2].type;
        cout << "Atom removed, -" << (rm_idx2 < rm_idx1 ? rm_idx2+1 : rm_idx2+2);
        cout << " " << cell_new.atm_list[rm_idx2].ele->sym_ << endl;
        cell_new.rm_atom(rm_idx2);
        return cell_new;
    }

    throw std::runtime_error("Neither add nor drop valid");
}


double DoubleMove::prefactor(const Cell &old, const Cell &trial) {
    double prefactor = phase_space_volume(old, trial);

    double add_weight_tot = 0.;
    for (auto elem : trial.ele_list)
        add_weight_tot += elem.p_add_;

    // Prob of choosing reverse move over forward move
    if (add_avail_) {
        // Combinatorial factor of two is not needed as it applies to both forward and reverse move
        // Forward move
        prefactor /= prob_add_
            * trial.ele_list[add_type_[0]].p_add_ / add_weight_tot
            * trial.ele_list[add_type_[1]].p_add_ / add_weight_tot;
        // Reverse move
        prefactor *= prob_drop_
            * trial.num_ele_each_remove[add_type_[0]]
            * (add_type_[0] != add_type_[1] ? trial.num_ele_each_remove[add_type_[1]]
                                            : trial.num_ele_each_remove[add_type_[1]]-1)
            / trial.num_atm_remove / (trial.num_atm_remove-1);
    } else if (drop_avail_) {
        // Forward move
        prefactor /= prob_drop_
            * old.num_ele_each_remove[drop_type_[0]]
            * (drop_type_[0] != drop_type_[1] ? old.num_ele_each_remove[drop_type_[1]]
                                              : old.num_ele_each_remove[drop_type_[1]]-1)
            / old.num_atm_remove / (old.num_atm_remove-1);
        // Reverse move
        prefactor *= prob_add_
            * old.ele_list[drop_type_[0]].p_add_ / add_weight_tot
            * old.ele_list[drop_type_[1]].p_add_ / add_weight_tot;
    } else {
        throw std::runtime_error("Previous move not recorded");
    }

    return prefactor;
}


// ----------------------------------------------------------------
// mc implementation
// ----------------------------------------------------------------

static std::map<std::string, double> read_action_block(std::istream &in) {
    const string label_act_p_start = "begin_action_probability";
    const string label_act_p_end = "end_action_probability";
    std::map<std::string, double> action_weight;

    // Find the block
    std::string line;
    in.seekg(ios::beg);
    in.clear();
    while (getline(in, line)) {
        if (line.find(label_act_p_start) != string::npos)
            break;
    }
    if (in.eof())
        throw std::runtime_error("Could not find action probabilities");

    while (getline(in, line)) {
        if (line.find(label_act_p_end) != std::string::npos)
            break;

        std::stringstream ss;
        ss << line;
        std::string label;
        ss >> label;
        // Lower case the string (note, crashes on non-ascii encodings)
        std::transform(label.begin(), label.end(), label.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        ss >> action_weight[label];
    }

    if (in.eof())
        throw std::runtime_error("Did not find end of action_probability block");

    return action_weight;
}


static std::vector<std::shared_ptr<MCMove>> init_actions_from_weights(
    const std::map<std::string, double> &action_weight
) {
    if (action_weight.at("add") < EPS || action_weight.at("drop") < EPS) {
        throw std::runtime_error("Add and Drop moves must both have positive weight");
    }
    auto ad_mov = std::make_shared<AddDropMove>(action_weight.at("add"), action_weight.at("drop"));
    std::vector<std::shared_ptr<MCMove>> moves{ad_mov};

    if (action_weight.find("swap") != action_weight.end()
        && action_weight.at("swap") > EPS
    ) {
        auto sp_mov = std::make_shared<SwapMove>(action_weight.at("swap"));
        moves.push_back(std::move(sp_mov));
    }

    if (action_weight.find("double") != action_weight.end()
        && action_weight.at("double") > EPS
    ) {
        auto db_mov = std::make_shared<DoubleMove>(action_weight.at("double"));
        moves.push_back(std::move(db_mov));
    }

    return moves;
}


void mc :: read_from_in(std::istream& in) {
    const string label_act_p = "action_probability";
    const string label_act_p_start = "begin_action_probability";

    read(in,"max_iter",'=',max_iter);
    read(in,"temperature",'=',temperature);
    read(in,"if_test",'=',if_test);

    // get action probability
    std::map<std::string, double> action_weight;
    string tmp;
    in.seekg(ios::beg);
    in.clear();
    while (getline(in, tmp)) {
        if (tmp.find(label_act_p_start) != string::npos) {
            action_weight = read_action_block(in);
            break;
        }
        if (tmp.find(label_act_p) != string::npos) {
            for (auto k : {"add", "drop", "swap"}) {
                in >> action_weight[k];
            }
            break;
        }
    }
    if (in.eof()) {
        throw std::runtime_error("Could not find action probabilities");
    }
    in.clear();
    in.seekg(ios::beg);

    // Initialize move objects
    moves_ = init_actions_from_weights(action_weight);

    // initialize other parameters
    act_type_ = -1;
    opt_e = 0;
}


std::shared_ptr<MCMove> mc::choose_next_move(const Cell &cell) {
    double action_weight[moves_.size()];

    //=======================================
    // check if each action is accessable, and adjust p if not
    cout << "Begin adjust weight of actions:" << endl;
    cout << "    New weight of each action is:" << endl;
    //---------------------------------------

    for (int i = 0; i < moves_.size(); i++) {
        auto move = moves_[i];
        if (move->available(cell)) {
            action_weight[i] = move->weight_;
        } else {
            action_weight[i] = 0.;
        }
        cout << "   " << move->name_ << ": " << action_weight[i];
    }

    cout << endl << "End adjust weight of actions" << endl << endl;

    //=======================================
    // decide which action to choose
    double total_weight = 0.;
    for (auto wt : action_weight) {
        total_weight += wt;
    }

    // exit if no actions are allowed
    if (fabs(total_weight) <= 1e-10) {
        throw std::runtime_error("No actions are available");
    }

    // Choose the action
    total_weight *= gcmc::rand_uniform();
    for (int i = 0; i < moves_.size(); i++) {
        total_weight -= action_weight[i];
        if (total_weight < 0) {
            act_type_ = i;
            return moves_[i];
        }
    }

    throw std::runtime_error("No action chosen");
}


std::shared_ptr<MCMove> mc::get_last_move() const {
    return moves_[act_type_];
}


void mc::create_new_structure(const Cell c_old, Cell& c_new) {
    // Choose move
    std::shared_ptr<MCMove> move = choose_next_move(c_old);

    // start applying change
    c_new = move->get_new_structure(c_old);
}

void mc :: save_opt_structure(const Cell c_new) {
    opt_e = c_new.energy;
    for(int t1=0; t1<c_new.num_atm; t1++)
        opt_e -= c_new.atm_list[t1].ele->mu_;
    opt_c = c_new;
    cout<<"Initialized the minimum seeker to the starting structure"<<endl;
}


static double formation_energy(const Cell &cell) {
    double energy = cell.energy;
    for (int i = 0; i < cell.num_atm; i++) {
        energy -= cell.atm_list[i].ele->mu_;
    }
    return energy;
}


bool mc :: check_if_accept(Cell& c_old, Cell& c_new) {
    // calculate formation energy
    e_old_ = formation_energy(c_old);
    e_trial_ = formation_energy(c_new);

    // Unclear what this does
    c_old.update_temperature(temperature);
    c_new.update_temperature(temperature);

    //==============================================
    cout<<endl<<"Evaluating whether to accept new structure"<<endl;
    if (fabs(c_new.energy) < 1e-9) {
        accept_ = false;
        cout << "    Warning: SCF of new structure does not converge" << endl;
        cout << "    ====Rejected====" << endl;
        return 0;
    }

    // Check optimal structure
    if (e_trial_ < opt_e) {
        opt_e = e_trial_;
        opt_c = c_new;
        cout<<"    New structure has by far the lowest formation energy, best structure updated"<<endl;
    }

    // Ratio of probability densities (independent of move type)
    double exp_main = exp(-(e_trial_-e_old_)/temperature/KB_EV);

    // Ratio of proposal densities
    double exp_pre = get_last_move()->prefactor(c_old, c_new);

    // start evaluating whether to accept or reject
    double accept_prob = exp_pre * exp_main;
    cout << "    Pre. factor: " << exp_pre;
    cout << "    Exp. factor: " << exp_main;
    cout << "    Accept prob: " << accept_prob << endl;

    if (1 < accept_prob) {
        accept_ = true;
        cout<<"    ====Accepted===="<<endl;
        cout<<"    New structure has lower volume-factored-in formation energy"<<endl;
    } else if (gcmc::rand_uniform() < accept_prob) {
        accept_ = true;
        cout<<"    ----Accepted----"<<endl;
        cout<<"    New structure has highter volume-factored-in formation energy"<<endl;
    } else {
        accept_ = false;
        cout<<"    ----Rejected----"<<endl;
    }

    cout<<"Best formation energy is "<<setw(20)<<setprecision(9)<<opt_e<<" eV"<<endl;
    return accept_;


}


void mc::log_header(std::ostream &log, const Cell &cell) const {
    log << setw(9) << "Iteration";
    for (int t1 = 0; t1 < cell.num_ele; t1++)
        log << setw(4) << cell.ele_list[t1].sym_;
    log << setw(20) << "DFT_E";
    log << setw(20) << "Trial_E_f";
    log << setw(20) << "Previous_E_f";
    log << setw(20) << "Accept_E_f";
    log << setw(20) << "Optimal_E_f";
    log << setw(14) << "If_accept" << endl;
}


void mc::log_step(std::ostream &log, const Cell &cell, int iter) const {
    log << setw(9) << iter;
    for(int t1 = 0; t1 < cell.num_ele; t1++)
        log << setw(4) << cell.num_ele_each[t1];
    log << fixed << setprecision(9) << setw(20) << cell.energy;
    log << setw(20) << e_trial_;
    log << setw(20) << e_old_;
    if (accept_) {
        log << setw(20) << e_trial_;
    } else {
        log << setw(20) << e_old_;
    }

    log << setw(20) << opt_e;
    log << setw(14) << accept_ << endl;
}


void mc :: print() {
    cout<<"Max iteration: "<<max_iter<<endl;
    cout<<"Simulation temperature: "<<temperature<<endl;
    cout<<"If running test: "<<if_test<<endl;
    cout<<"Action probability: "<<endl;
    for (auto mov : moves_)
        cout << "\t" << mov->name_ << ": " << mov->weight_ << endl;
    cout << endl;
}
