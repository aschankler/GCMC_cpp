#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <random>
#include "mc.h"
#include "auxiliary.h"
#include "rng.h"


// In eV/K
#define KB_EV 8.6173303e-5

using namespace std;

void mc :: read_from_in(ifstream& in) {
    string label_act_p = "begin_action_probability";
    string tmp;
    stringstream ss;

    read(in,"max_iter",'=',max_iter);
    read(in,"temperature",'=',temperature);
    read(in,"if_test",'=',if_test);

    // get action probability
    in.seekg(ios::beg);
    in.clear();
    while(getline(in, tmp))
        if(tmp.find(label_act_p) != string::npos)
            break;
    for(int t1=0; t1<num_act; t1++)
        in>>act_p[t1];
    in.clear();
    in.seekg(ios::beg);

    // initialize other parameters
    act_type = -1;
    opt_e = 0;
}


int choose_add(const Cell &cell_old, vec &pos_add, int &ele_type) {
    const int num_trial = 1000;

    // Choose species to add
    double add_weight_sum = 0;
    for (int t1 = 0; t1 < cell_old.num_ele; t1++)
        add_weight_sum += cell_old.ele_list[t1].p_add;

    // Choose random number from [0, add_weight_sum)
    ele_type = cell_old.num_ele;
    std::uniform_real_distribution<double> ele_dist(0., add_weight_sum);
    add_weight_sum = ele_dist(gcmc::rng);
    for (int t1 = 0; t1 < cell_old.num_ele; t1++) {
        add_weight_sum -= cell_old.ele_list[t1].p_add;
        if(add_weight_sum <= 0) {
            ele_type = t1;
            break;
        }
    }

    element ele_tmp = cell_old.ele_list[ele_type];

    // Choose position to add
    for (int t1 = 0; t1 < num_trial; t1++) {
        // Generate new position
        pos_add.rand();  // add in [0, 1)^3; crystal coordinates

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
        if (pos_add.x[2] < cell_old.h_min || pos_add.x[2] > cell_old.h_max)
            continue;

        // Filter coordination rule
        if (cell_old.if_vc_relax) {
            double r_tmp = cell_old.min_distance(pos_add);
            if (r_tmp > ele_tmp.r_min && r_tmp < ele_tmp.r_max) {
                // Accept this position
                return 1;
            }
        }
    }

    // Did not find a suitable posititon
    return 0;
}

int choose_drop(const Cell &cell_old, int &rm_idx) {
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


int choose_swap(const Cell &cell_old, int &idx1, int &idx2) {
    std::uniform_int_distribution<> dist(0, cell_old.num_atm_move-1);

    // find the first atom to switch
    int atm_id_tmp = dist(gcmc::rng);
    for (int t1 = 0; t1 < cell_old.num_atm; t1++) {
        if (cell_old.atm_list[t1].if_move >= 1) {
            // Atom is movable
            if (atm_id_tmp == 0) {
                idx1 = t1;
                break;
            }
            atm_id_tmp--;
        }
    }

    // find the second atom to switch
    int num_move_distinct = cell_old.num_atm_move - cell_old.num_ele_each_move[cell_old.atm_list[idx1].type];
    atm_id_tmp = dist(gcmc::rng, decltype(dist)::param_type(0, num_move_distinct-1));
    for (int t1 = 0; t1 < cell_old.num_atm; t1++) {
        if (cell_old.atm_list[t1].if_move >= 1 && cell_old.atm_list[t1].type != cell_old.atm_list[idx1].type) {
            // Atom is movable
            if (atm_id_tmp == 0) {
                idx2 = t1;
                return 1;
            }
            atm_id_tmp--;
        }
    }
    return 0;
}


void mc :: create_new_structure(const Cell c_old, Cell& c_new) {
    c_new = c_old;
    double tmp_p[num_act];
    // for examine add
    vec add_pos;
    int ele_type_add;

    // for exam remove and swap
    int rm_idx;
    int swp_idx1, swp_idx2;
    int num_movable_ele;

    double act_p_sum;

    num_atm_each_change_.resize(c_new.num_ele);
    for (int t1 = 0; t1 < c_new.num_ele; t1++)
        num_atm_each_change_[t1] = 0;

    //=======================================
    // check if each action is accessable, and adjust p if not
    cout<<"Begin adjust weight of actions:"<<endl;
    //---------------------------------------
    // check for add
    if (choose_add(c_old, add_pos, ele_type_add)) {
        tmp_p[0] = act_p[0];
    } else {
        cout << "    Can not find site to add atom, weight of choosing to add is set to 0" << endl;
        tmp_p[0] = 0.;
    }

    //---------------------------------------
    // check for remove
    if (c_new.num_atm_remove > 0) {
        tmp_p[1] = act_p[1];
    } else {
        tmp_p[1] = 0;
        cout << "    Can not find removable elements, weight of choosing to remove set to 0" << endl;
    }

    // check for swap
    num_movable_ele = 0;
    for (int t1 = 0; t1 < c_new.num_ele; t1++) {
        if (c_new.num_ele_each_move[t1] > 0)
            num_movable_ele++;
    }
    if (num_movable_ele >= 2) {
        tmp_p[2] = act_p[2];
    } else {
        tmp_p[2] = 0;
        cout<<"    Can not find more than one movable elements, weight of choosing to swap is set to 0"<<endl;
    }

    cout << "    New weight of each action is:" << endl;
    cout << "    Add: " << tmp_p[0];
    cout << "    Remove: " << tmp_p[1];
    cout << "    Swap: " << tmp_p[2] << endl;
    cout << "End adjust weight of actions" << endl << endl;

    //=======================================
    // create new structure
    // decide which action to choose
    act_p_sum = 0;
    for(int t1 = 0; t1 < num_act; t1++)
        act_p_sum += tmp_p[t1];

    // exit if no actions are allowed
    if (fabs(act_p_sum) <= 1e-10) {
        cout<<"Error: No actions are legal"<<endl;
        exit(EXIT_FAILURE);
    }

    act_p_sum = (double)rand()/RAND_MAX * act_p_sum;
    act_type_ = -1;
    for (int t1 = 0; t1 < num_act; t1++) {
        act_p_sum -= tmp_p[t1];
        if (act_p_sum < 0) {
            act_type_ = t1;
            break;
        }
    }


    // start applying change
    //=======================================
    switch(act_type_) {
    //---------------------------------------
    // add
    case 0: {
        num_atm_each_change_[ele_type_add] = 1;
        c_new.ad_atom(add_pos, ele_type_add);
        cout << "Atom added, +"<<c_new.num_atm<<" "<<c_new.ele_list[ele_type_add].sym<<", position is "<<add_pos<<endl;
        break;
    }
    //---------------------------------------
    // remove
    case 1: {
        choose_drop(c_new, rm_idx);
        num_atm_each_change_[c_old.atm_list[rm_idx].type] = -1;
        c_new.rm_atom(rm_idx);
        cout<<"Atom removed, -"<<rm_idx+1<<" "<<c_old.atm_list[rm_idx].ele->sym<<endl;
        break;
    }
    //---------------------------------------
    // swap
    case 2: {
        //int iter_swap = rand()%c_old.num_atm+1;
        int iter_swap = 1; // only swap one pair
        cout<<"Perform swap "<<iter_swap<<" times:"<<endl;
        for(int t_swap=0; t_swap<iter_swap; t_swap++) {
            choose_swap(c_new, swp_idx1, swp_idx2);
            c_new.sp_atom(swp_idx1, swp_idx2);
            cout<<"Atoms swapped, "<<swp_idx1+1<<" "<<c_old.atm_list[swp_idx1].ele->sym<<" <--> "<<swp_idx2+1<<" "<<c_old.atm_list[swp_idx2].ele->sym<<endl;
        }
        break;
    }
    //---------------------------------------
    default: {
        cout<<"Error: Undifined action type: "<<act_type<<endl;
        exit(EXIT_FAILURE);
    }
    }
}

void mc :: save_opt_structure(const Cell c_new) {
    opt_e = c_new.energy;
    for(int t1=0; t1<c_new.num_atm; t1++)
        opt_e -= c_new.atm_list[t1].ele->mu;
    opt_c = c_new;
    cout<<"Initialized the minimum seeker to the starting structure"<<endl;
}

int mc :: check_if_accept(Cell& c_old, Cell& c_new) {
    double exp_pre, exp_main;

    // calculate formation energy
    e_old_ = c_old.energy;
    e_trial_ = c_new.energy;
    for (int t1 = 0; t1 < c_old.num_atm; t1++)
        e_old_ -= c_old.atm_list[t1].ele->mu;
    for (int t1 = 0; t1 < c_new.num_atm; t1++)
        e_trial_ -= c_new.atm_list[t1].ele->mu;

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

    switch (act_type_) {
    //---------------------------------------
    // add,remove, and swap
    case 0:
    case 1:
    case 2: {
        // calculate prefactor and exp
        exp_pre = 1;
        for (int t1 = 0; t1 < c_old.num_ele; t1++) {
            c_old.ele_list[t1].update_tb(temperature);
            c_new.ele_list[t1].update_tb(temperature);
            exp_pre *= pow(c_new.get_volume(), num_atm_each_change_[t1]);
            exp_pre *= factor(c_old.num_ele_each[t1]) / factor(num_atm_each_change_[t1] + c_old.num_ele_each[t1]);
            exp_pre /= pow(c_old.ele_list[t1].tb, 3*num_atm_each_change_[t1]);
        }

        // if allow V to change, multiply exp_pre by (V_new/V_old)^N_tot(old)
        if(c_new.if_change_v)
            exp_pre *= pow((c_new.get_volume()/c_old.get_volume()), c_old.num_atm);

        exp_main = exp(-(e_trial_-e_old_)/temperature/KB_EV);

        cout<<"    Pre. factor: "<<exp_pre<<"    exp. factor: "<<exp_main<<"    total factor: "<<exp_pre*exp_main<<endl;

        // start evaluating whether to accept or reject
        if (1 < exp_pre*exp_main) {
            accept_ = true;
            cout<<"    ====Accepted===="<<endl;
            cout<<"    New structure has lower volume-factored-in formation energy"<<endl;
        } else if ((double)rand()/RAND_MAX < exp_pre*exp_main) {
            accept_ = true;
            cout<<"    ----Accepted----"<<endl;
            cout<<"    New structure has highter volume-factored-in formation energy"<<endl;
        } else {
            accept_ = false;
            cout<<"    ----Rejected----"<<endl;
        }
        break;
    }
    //---------------------------------------
    default: {
        cout<<"Error: Invalid action number: "<<act_type<<endl;
        exit(EXIT_FAILURE);
    }
    }

    cout<<"Best formation energy is "<<setw(20)<<setprecision(9)<<opt_e<<" eV"<<endl;
    return accept_;
}

int mc :: factor(int n) {
    if (n < 0) {
        cout<<"Error: Number of atoms should not be negative"<<endl;
        exit(EXIT_FAILURE);
    }
    if (n==0)
        return 1;
    else
        return factor(n-1);
}


void mc::log_header(std::ostream &log, const Cell &cell) const {
    log << setw(9) << "Iteration";
    for (int t1 = 0; t1 < cell.num_ele; t1++)
        log << setw(4) << cell.ele_list[t1].sym;
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
    for(int t1=0; t1<num_act; t1++)
        cout<<'\t'<<act_p[t1];
    cout<<endl;
}
