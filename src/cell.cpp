#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "atom.h"
#include "cell.h"
#include "vec.h"
#include "auxiliary.h"

#define  RY_EV  13.605693009

using namespace std;


CellControlParams cell_control_from_in(std::istream &in) {
    double hmin, hmax;

    read(in, "h_min", '=', hmin);
    read(in, "h_max", '=', hmax);

    CellControlParams params{hmin, hmax};

    read_opt(in, "a_min", '=', params.a_min_);
    read_opt(in, "a_max", '=', params.a_max_);
    read_opt(in, "b_min", '=', params.b_min_);
    read_opt(in, "b_max", '=', params.b_max_);

    read_opt(in, "if_vc_relax", '=', params.if_vc_relax_);
    read_opt(in, "if_change_v", '=', params.if_change_v_);

    return params;
}


Cell cell_from_in(istream& in) {
    string label_ele = "begin_elements";
    string label_lat = "begin_lattice";
    string label_atm = "begin_atom_positions";
    string tmp;
    stringstream ss;

    Cell cell{cell_control_from_in(in)};

    read(in, "num_ele", '=', cell.num_ele);
    read(in, "num_atm", '=', cell.num_atm);

    // generate element list
    while(getline(in,tmp))
        if(tmp.find(label_ele) != string::npos)
            break;
    for (int t1 = 0; t1 < cell.num_ele; t1 ++) {
        getline(in, tmp);
        cell.ele_list.push_back(element_from_input(tmp));
    }
    in.clear();
    in.seekg(ios::beg);

    // generate atom list
    while(getline(in,tmp))
        if(tmp.find(label_atm) != string::npos)
            break;
    for(int t1=0; t1<cell.num_atm; t1++) {
        getline(in,tmp);
        cell.atm_list.push_back(atom_from_input(tmp, cell.ele_list));
    }
    in.clear();
    in.seekg(ios::beg);

    // get lattice parameter
    while(getline(in,tmp))
        if(tmp.find(label_lat) != string::npos)
            break;
    in>>cell.latt[0]>>cell.latt[1]>>cell.latt[2];
    in.clear();
    in.seekg(ios::beg);

    cell.update();

    return cell;
}

void Cell::read_output(int calculator_type) {
    switch (calculator_type) {
    case 1: { //QE
        read_from_qe(*this);
        break;
    }
    case 2: { //vasp
        read_from_vasp(*this);
        break;
    }
    default: {
        cout<<"Error: Invalid calculator type!"<<endl;
        exit(EXIT_FAILURE);
    }
    }

    count_move_atoms();
    update_volume();
    update_lat_inv();
}


void Cell::zero_force() {
    for(int t1 = 0; t1 < num_atm; t1++)
        atm_list[t1].force_ = atm_list[t1].force_ * 0;
}


void read_from_qe(Cell &cell_out) {
    ifstream in("qe.out");

    string label_num_atm = "number of atoms/cell";
    string label_force = "Forces acting on atoms";
    string label_cell = "CELL_PARAMETERS";
    string label_position = "ATOMIC_POSITIONS";
    string label_energy = "!    total energy";
    string label_final_energy = "Final energy";
    string label_final_enthalpy = "Final enthalpy";
    const string label_no_scf_conv = "convergence NOT achieved";

    //string label_final_position = "Begin final coordinates";
    string tmp;
    stringstream ss;
    int num_tmp;
    int pos_tmp, energy_pos;
    bool no_scf_conv = false;

    // check number of atoms matches
    read(in, label_num_atm, '=', num_tmp);
    if(num_tmp != cell_out.num_atm) {
        cout<<"Error: Number of atoms in qe.out does not match record"<<endl;
        exit(EXIT_FAILURE);
    }

    // find number of iteration
    for (num_tmp = 0; !in.eof(); getline(in, tmp)) {
        if (tmp.find(label_position) != string::npos)
            num_tmp++;
        if (tmp.find(label_no_scf_conv) != string::npos)
            no_scf_conv = true;
    }
    in.clear();
    in.seekg(ios::beg);

    // check if SCF converges
    if (num_tmp == 0 or no_scf_conv) {
        cout<<"Warning: SCF does not converge, set energy and force to 0"<<endl;
        cell_out.energy = 0;
        cell_out.zero_force();
        in.clear();
        in.seekg(ios::beg);
        return;
    }

    // save energy, forces, (cell parameters), and positions
    for(int t1 = 0; t1 < num_tmp;) {
        pos_tmp = in.tellg();
        getline(in, tmp);
        if (tmp.find(label_energy) != string::npos) {
            energy_pos = pos_tmp;
            t1++;
        }
    }

    // energy
    ss << (tmp);
    getline(ss, tmp, '=');
    ss >> cell_out.energy;
    cell_out.energy *= RY_EV;
    ss.str("");
    ss.clear();

    // forces
    while (!in.eof() and tmp.find(label_force) == string::npos)
        getline(in, tmp);
    if (in.eof()) {
        cout << "Could not find forces" << endl;
        cell_out.zero_force();
    } else {
        getline(in, tmp);
        for(int t1 = 0; t1 < cell_out.num_atm; t1++) {
            getline(in, tmp, '=');
            in >> cell_out.atm_list[t1].force_;
            getline(in, tmp);
        }
    }

    // cell parameters
    if (cell_out.control.if_vc_relax_) {
        while(tmp.find(label_cell) == string::npos)
            getline(in, tmp);
        for(int t1 = 0; t1 < 3; t1++) {
            in >> cell_out.latt[t1];
            getline(in, tmp);
        }
    }

    // position
    in.clear();
    in.seekg(energy_pos);
    while(!in.eof() and tmp.find(label_position) == string::npos)
        getline(in, tmp);
    if (in.eof()) {
        cout << "Error: could not find positions" << endl;
        exit(EXIT_FAILURE);
    } else {
        for(int t1 = 0; t1 < cell_out.num_atm; t1++) {
            in >> tmp >> cell_out.atm_list[t1].pos_;
            getline(in, tmp);
        }
    }

    // check if relax converged
    in.clear();
    in.seekg(ios::beg);
    while(tmp.find(label_final_energy) == string::npos && tmp.find(label_final_enthalpy) == string::npos && !in.eof())
        getline(in, tmp);
    if(in.eof())
        cout<<"Warning: Structure not fully relaxed"<<endl;
    in.clear();
    in.seekg(ios::beg);

    in.close();
}

void read_from_vasp(Cell &cell_out) {
    ifstream in("vasp.out");
    ifstream coord_in("CONTCAR");
    string label_energy = "F=";
    string label_relaxed = "reached required accuracy";
    string label_atm = "Direct";
    stringstream ss;
    string tmp;
    int num_tmp;

    // save energy
    num_tmp = 0;
    while(getline(in, tmp)) {
        if (tmp.find(label_energy) != string::npos) {
            ss << (tmp);
            getline(ss, tmp, '=');
            ss >> cell_out.energy;
            ss.str("");
            ss.clear();
            num_tmp++;
        } else if (tmp.find(label_relaxed) != string::npos)
            break;
    }

    // check if SCF converges
    if (num_tmp == 0) {
        cout<<"Warning: SCF does not converge, set energy and force to 0"<<endl;
        cell_out.energy = 0;
        cell_out.zero_force();
        return;
    }

    // check if relaxed
    if (in.eof())
        cout<<"Warning: Structure not fully relaxed"<<endl;
    in.close();

    // save cell parameters
    getline(coord_in, tmp);
    getline(coord_in, tmp);
    coord_in >> cell_out.latt[0];
    coord_in >> cell_out.latt[1];
    coord_in >> cell_out.latt[2];

    // structure file check
    getline(coord_in, tmp);
    getline(coord_in, tmp);
    for(int t1 = 0; t1 < cell_out.num_ele; t1++) {
        coord_in >> num_tmp;
        if(num_tmp != cell_out.num_ele_each[t1]) {
            throw std::runtime_error("Wrong number of atoms in CONTCAR!");
        }
    }

    // save atomic positions
    while(getline(coord_in, tmp))
        if(tmp.find(label_atm) != string::npos)
            break;

    if(coord_in.eof()) {
        cout<<"Warning: Final structure not found, set energy and force to 0"<<endl;
        cell_out.energy = 0;
        cell_out.zero_force();
        return;
    } else {
        vec frac_coord, atomic_coord;
        for (int i_ele = 0; i_ele < cell_out.num_ele; i_ele++) {
            for(auto& m2:cell_out.atm_list) {
                if (i_ele == m2.type_) {
                    coord_in >> frac_coord;
                    getline(coord_in, tmp);
                    atomic_coord = cell_out.latt[0]*frac_coord.x[0] + cell_out.latt[1]*frac_coord.x[1] + cell_out.latt[2]*frac_coord.x[2];
                    m2.pos_ = atomic_coord;
                }
            }
        }
        cout<<"Warning: Forces not yet implemented for VASP, set forces to 0"<<endl;
        cell_out.zero_force();
    }
}

void Cell::write_axsf(ofstream& out) const {
    out<<"ANIMSTEPS 10000"<<endl;
    out<<"CRYSTAL"<<endl;
    out<<"PRIMVEC"<<endl;
    out<<latt[0]<<endl;
    out<<latt[1]<<endl;
    out<<latt[2]<<endl;
}

void Cell::write_axsf(ofstream& out,int iter) const {
    out<<"PRIMCOORD "<<iter<<endl;
    out<<num_atm<<" 1"<<endl;
    for(int t1=0; t1<num_atm; t1++) {
        out << setw(2) << atom_type(t1).sym_;
        out<<atm_list[t1].pos_<<atm_list[t1].force_<<endl;
    }
}

void Cell::write_xsf(ofstream& out,int iter) const {
    out<<"PRIMVEC"<<endl;
    out<<latt[0]<<endl;
    out<<latt[1]<<endl;
    out<<latt[2]<<endl;
    out<<"PRIMCOORD "<<iter<<endl;
    out<<num_atm<<" 1"<<endl;
    for(int t1=0; t1<num_atm; t1++) {
        out << setw(2) << atom_type(t1).sym_;
        out<<atm_list[t1].pos_<<atm_list[t1].force_<<endl;
    }
}


void Cell::update() {
    update_volume();
    update_lat_inv();
    count_move_atoms();
}


void Cell::count_move_atoms() {
    num_ele_each.resize(num_ele);
    num_ele_each_move.resize(num_ele);
    num_ele_each_remove.resize(num_ele);
    num_atm_move = 0;
    num_atm_remove = 0;
    for(int t1=0; t1<num_ele; t1++) {
        num_ele_each[t1] = 0;
        num_ele_each_move[t1] = 0;
        num_ele_each_remove[t1] = 0;
    }
    // start counting
    for(int t1=0; t1<num_atm; t1++) {
        num_ele_each[atm_list[t1].type_]++;
        switch(atm_list[t1].if_move_) {
        // not movable
        case 0: {
            break;
        }
        // movable not removable
        case 1: {
            num_atm_move++;
            num_ele_each_move[atm_list[t1].type_]++;
            break;
        }
        // all free
        case 2: {
            num_atm_move++;
            num_atm_remove++;
            num_ele_each_move[atm_list[t1].type_]++;
            num_ele_each_remove[atm_list[t1].type_]++;
            break;
        }
        default: {
            cout<<"Error: Atom "<<t1+1<<' '<<atom_type(t1).sym_<<" does not have a vaild (re)movable flag"<<endl;
            exit(EXIT_FAILURE);
        }
        }
    }
}


void Cell::update_volume() {
    double h[3]= {0,0,control.h_max_};
    vec hh;

    hh = &h[0];
    if (control.if_vc_relax_)
        vol_ = (latt[0]^latt[1])*latt[2];
    else
        vol_ = (latt[0]^latt[1])*hh;
}


void Cell::update_lat_inv() {
    latt_inv[0] = (latt[1]^latt[2])/((latt[0]^latt[1])*latt[2]);
    latt_inv[1] = (latt[2]^latt[0])/((latt[0]^latt[1])*latt[2]);
    latt_inv[2] = (latt[0]^latt[1])/((latt[0]^latt[1])*latt[2]);
}


void Cell::update_temperature(double temperature) {
    for (int i = 0; i < num_ele; i++) {
        ele_list[i].update_tb(temperature);
    }
}


const vec Cell::to_crystal(const vec& pos) const {
    vec crystal_coord;
    crystal_coord.clean();

    // Note: latt_inv = (latt^-1)^T
    // crystal = latt_inv @ cart
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            crystal_coord[i] += latt_inv[i][j] * pos[j];
        }
    }
    return crystal_coord;
}


const vec Cell::from_crystal(const vec& pos) const {
    vec cart_coord;
    cart_coord.clean();

    // cart = crystal @ latt
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cart_coord[i] += pos[j] * latt[j][i];
        }
    }
    return cart_coord;
}


void Cell::ad_atom(vec pos, int ele_type) {
    Atom tmp;
    tmp.type_ = ele_type;
    tmp.pos_ = pos;
    tmp.force_ = pos*0;
    tmp.if_move_ = 2;
    atm_list.push_back(tmp);
    num_atm++;
    count_move_atoms();
}

void Cell::rm_atom(int ind_atm) {
    if (atm_list[ind_atm].if_move_ <= 1) {
        cout<<"Error: Can not remove atom "<<ind_atm+1<<' '<<atom_type(ind_atm).sym_<<", not removable"<<endl;
        exit(EXIT_FAILURE);
    }
    atm_list.erase(atm_list.begin() + ind_atm);
    num_atm--;
    count_move_atoms();
}

void Cell::sp_atom(int s1, int s2) {
    vec tmp;
    tmp = atm_list[s1].pos_;
    atm_list[s1].pos_ = atm_list[s2].pos_;
    atm_list[s2].pos_ = tmp;

    tmp = atm_list[s1].force_;
    atm_list[s1].force_ = atm_list[s2].force_;
    atm_list[s2].force_ = tmp;
}

void Cell::update_tb(double T) {
    for(int t1=0; t1<num_ele; t1++)
        ele_list[t1].update_tb(T);
}


/* Find shortest distance between `pos` and another atom in the cell */
double Cell::min_distance(vec pos, int& min_idx) const {
    double r_min, r_tmp;
    r_min = 1e10;

    for (int t1 = -1; t1 < 2; t1++) {
        for (int t2 = -1; t2 < 2; t2++) {
            for (int t3 = -1; t3 < 2; t3++) {
                for (int at_idx = 0; at_idx < num_atm; at_idx++) {
                    r_tmp = (
                        pos + latt[0]*t1 + latt[1]*t2 + latt[2]*t3 - atm_list[at_idx].pos_
                    ).norm();

                    if (r_min > r_tmp) {
                        r_min = r_tmp;
                        min_idx = at_idx;
                    }
                }
            }
        }
    }

    return r_min;
}

double Cell::min_distance(vec pos) const {
    int min_idx;
    return min_distance(pos, min_idx);
}

void Cell :: print() const {
    cout<<"Number of elements: "<<num_ele<<endl;
    cout<<"Number of atoms: "<<num_atm<<endl;
    cout<<"Threshold of height: "<<control.h_min_ << ", " << control.h_max_ << endl;
    cout<<"Energy: "<<energy<<endl;
    cout<<"Volume: " << get_volume() << endl;
    cout<<"Total removable atoms: "<<num_atm_remove<<endl;
    cout<<"Total movable atoms: "<<num_atm_move<<endl;
    cout<<"Atoms per elements: "<<endl;
    for(int t1=0; t1<num_ele; t1++)
        cout<<ele_list[t1].sym_<<'\t'<<num_ele_each[t1]<<endl;
    cout<<"Reovable atoms per elements: "<<endl;
    for(int t1=0; t1<num_ele; t1++)
        cout<<ele_list[t1].sym_<<'\t'<<num_ele_each_remove[t1]<<endl;
    cout<<"Movable atoms per elements: "<<endl;
    for(int t1=0; t1<num_ele; t1++)
        cout<<ele_list[t1].sym_<<'\t'<<num_ele_each_move[t1]<<endl;
//	cout<<"List of elements:"<<endl;
//	for(int t1=0; t1<num_ele; t1++)
//		ele_list[t1].print();
    cout<<"List of atoms:"<<endl;
    for(int t1=0; t1<num_atm; t1++)
        atm_list[t1].print();
}
