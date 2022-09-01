#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "atom.h"
#include "cell.h"
#include "vec.h"
#include "auxiliary.h"

using namespace std;

void cell :: read_from_in(ifstream& in) {
    string label_ele = "begin_elements";
    string label_lat = "begin_lattice";
    string label_atm = "begin_atom_positions";
    string tmp;
    stringstream ss;

    read(in,"num_ele",'=',num_ele);
    read(in,"num_atm",'=',num_atm);
    read(in,"h_min",'=',h_min);
    read(in,"h_max",'=',h_max);
    read_opt(in,"a_min",'=',a_min);
    read_opt(in,"a_max",'=',a_max);
    read_opt(in,"b_min",'=',b_min);
    read_opt(in,"b_max",'=',b_max);
    read(in,"if_vc_relax",'=',if_vc_relax);
    read(in,"if_change_v",'=',if_change_v);

    // generate element list
    ele_list.resize(num_ele);
    while(getline(in,tmp))
        if(tmp.find(label_ele) != string::npos)
            break;
    for(int t1=0; t1<num_ele; t1++) {
        getline(in,tmp);
        ele_list[t1].get_param(tmp);
    }
    in.clear();
    in.seekg(ios::beg);

    // generate atom list
    atm_list.resize(num_atm);
    while(getline(in,tmp))
        if(tmp.find(label_atm) != string::npos)
            break;
    for(int t1=0; t1<num_atm; t1++) {
        getline(in,tmp);
        atm_list[t1].line_from_in(tmp,ele_list);
    }
    in.clear();
    in.seekg(ios::beg);

    // get lattice parameter
    while(getline(in,tmp))
        if(tmp.find(label_lat) != string::npos)
            break;
    in>>latt[0]>>latt[1]>>latt[2];
    latt_inv[0] = (latt[1]^latt[2])/((latt[0]^latt[1])*latt[2]);
    latt_inv[1] = (latt[2]^latt[0])/((latt[0]^latt[1])*latt[2]);
    latt_inv[2] = (latt[0]^latt[1])/((latt[0]^latt[1])*latt[2]);
    in.clear();
    in.seekg(ios::beg);
}

void cell :: read_output(int calculator_type) {
    switch (calculator_type) {
    case 1: { //QE
        read_from_qe();
        break;
    }
    case 2: { //vasp
        read_from_vasp();
        break;
    }
    default: {
        cout<<"Error: Invalid calculator type!"<<endl;
        exit(EXIT_FAILURE);
    }
    }
}

void cell :: read_from_qe() {
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
    bool no_scf_conv = false;

    // check number of atoms matches
    read(in, label_num_atm, '=', num_tmp);
    if(num_tmp != num_atm) {
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
        energy = 0;
        for(int t1=0; t1<num_atm; t1++)
            atm_list[t1].force = atm_list[t1].force * 0;
        in.clear();
        in.seekg(ios::beg);
        return;
    }
    // save energy, forces, (cell parameters), and positions
    for(int t1=0; t1<num_tmp;) {
        getline(in,tmp);
        if (tmp.find(label_energy) != string::npos)
            t1++;
    }
    // energy
    ss<<(tmp);
    getline(ss,tmp,'=');
    ss>>energy;
    energy*=ry_ev;
    ss.str("");
    ss.clear();
    // forces
    while(tmp.find(label_force) == string::npos)
        getline(in,tmp);
    getline(in,tmp);
    for(int t1=0; t1<num_atm; t1++) {
        getline(in,tmp,'=');
        in>>atm_list[t1].force;
        getline(in,tmp);
    }
    // cell parameters
    if (if_vc_relax) {
        while(tmp.find(label_cell) == string::npos)
            getline(in,tmp);
        for(int t1=0; t1<3; t1++) {
            in>>latt[t1];
            getline(in,tmp);
        }
    }
    // position
    while(tmp.find(label_position) == string::npos)
        getline(in,tmp);
    for(int t1=0; t1<num_atm; t1++) {
        in>>tmp>>atm_list[t1].pos;
        getline(in,tmp);
    }
    // check if relax converged
    in.clear();
    in.seekg(ios::beg);
    while(tmp.find(label_final_energy) == string::npos && tmp.find(label_final_enthalpy) == string::npos && !in.eof())
        getline(in,tmp);
    if(in.eof())
        cout<<"Warning: Structure not fully relaxed"<<endl;
    in.clear();
    in.seekg(ios::beg);

    in.close();
}

void cell :: read_from_vasp() {
    ifstream in("vasp.out");
    ifstream coord_in("CONTCAR");
    string label_energy = "F=";
    string label_relaxed = "reached required accuracy";
    string label_atm = "Direct";
    stringstream ss;
    string tmp;
    int num_tmp;

    // save energy
    num_tmp=0;
    while(getline(in,tmp)) {
        if (tmp.find(label_energy) != string::npos) {
            ss << (tmp);
            getline(ss,tmp,'=');
            ss >> energy;
            ss.str("");
            ss.clear();
            num_tmp++;
        } else if (tmp.find(label_relaxed) != string::npos)
            break;
    }
    // check if SCF converges
    if (num_tmp == 0) {
        cout<<"Warning: SCF does not converge, set energy and force to 0"<<endl;
        energy = 0;
        for(int t1=0; t1<num_atm; t1++)
            atm_list[t1].force = atm_list[t1].force * 0;
        return;
    }
    // check if relaxed
    if (in.eof())
        cout<<"Warning: Structure not fully relaxed"<<endl;
    in.close();

    // save cell parameters
    getline(coord_in,tmp);
    getline(coord_in,tmp);
    coord_in>>latt[0];
    coord_in>>latt[1];
    coord_in>>latt[2];
    get_volume();
    // structure file check
    getline(coord_in,tmp);
    getline(coord_in,tmp);
    for(int t1=0; t1<num_ele; t1++) {
        coord_in>>num_tmp;
        if(num_tmp != num_ele_each[t1]) {
            cout<<"Error: Wrong number of atoms in CONTCAR!"<<endl;
            exit(EXIT_FAILURE);
        }
    }
    // save atomic positions
    while(getline(coord_in,tmp))
        if(tmp.find(label_atm) != string::npos)
            break;
    if(coord_in.eof()) {
        cout<<"Warning: Final structure not found, set energy and force to 0"<<endl;
        energy = 0;
        for(int t1=0; t1<num_atm; t1++)
            atm_list[t1].force = atm_list[t1].force * 0;
        return;
    } else {
        vec frac_coord, atomic_coord;
        for(auto m1:ele_list) {
            for(auto& m2:atm_list) {
                if(m1.sym == m2.ele->sym) {
                    coord_in>>frac_coord;
                    getline(coord_in,tmp);
                    atomic_coord = latt[0]*frac_coord.x[0]+latt[1]*frac_coord.x[1]+latt[2]*frac_coord.x[2];
                    m2.pos = atomic_coord;
                }
            }
        }
        cout<<"Warning: Forces not yet implemented for VASP, set forces to 0"<<endl;
        for(int t1=0; t1<num_atm; t1++)
            atm_list[t1].force = atm_list[t1].force * 0;
    }
}

void cell :: write_axsf(ofstream& out) {
    out<<"ANIMSTEPS 10000"<<endl;
    out<<"CRYSTAL"<<endl;
    out<<"PRIMVEC"<<endl;
    out<<latt[0]<<endl;
    out<<latt[1]<<endl;
    out<<latt[2]<<endl;
}

void cell :: write_axsf(ofstream& out,int iter) {
    out<<"PRIMCOORD "<<iter<<endl;
    out<<num_atm<<" 1"<<endl;
    for(int t1=0; t1<num_atm; t1++) {
        out<<setw(2)<<atm_list[t1].ele->sym;
        out<<atm_list[t1].pos<<atm_list[t1].force<<endl;
    }
}

void cell :: write_xsf(ofstream& out,int iter) {
    out<<"PRIMVEC"<<endl;
    out<<latt[0]<<endl;
    out<<latt[1]<<endl;
    out<<latt[2]<<endl;
    out<<"PRIMCOORD "<<iter<<endl;
    out<<num_atm<<" 1"<<endl;
    for(int t1=0; t1<num_atm; t1++) {
        out<<setw(2)<<atm_list[t1].ele->sym;
        out<<atm_list[t1].pos<<atm_list[t1].force<<endl;
    }
}

void cell :: count_move_atoms() {
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
        num_ele_each[atm_list[t1].type]++;
        switch(atm_list[t1].if_move) {
        // not movable
        case 0: {
            break;
        }
        // movable not removable
        case 1: {
            num_atm_move++;
            num_ele_each_move[atm_list[t1].type]++;
            break;
        }
        // all free
        case 2: {
            num_atm_move++;
            num_atm_remove++;
            num_ele_each_move[atm_list[t1].type]++;
            num_ele_each_remove[atm_list[t1].type]++;
            break;
        }
        default: {
            cout<<"Error: Atom "<<t1+1<<' '<<atm_list[t1].ele->sym<<" does not have a vaild (re)movable flag"<<endl;
            exit(EXIT_FAILURE);
        }
        }
    }
}

double cell :: get_volume() {
    double h[3]= {0,0,h_max};
    vec hh;
    hh=&h[0];
    if (if_vc_relax)
        vol = (latt[0]^latt[1])*latt[2];
    else
        vol = (latt[0]^latt[1])*hh;

    latt_inv[0] = (latt[1]^latt[2])/((latt[0]^latt[1])*latt[2]);
    latt_inv[1] = (latt[2]^latt[0])/((latt[0]^latt[1])*latt[2]);
    latt_inv[2] = (latt[0]^latt[1])/((latt[0]^latt[1])*latt[2]);
    return vol;
}


const vec cell::to_crystal(const vec& pos) const {
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


const vec cell::from_crystal(const vec& pos) const {
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


void cell :: ad_atom(vec pos, int ele_type) {
    atom tmp;
    tmp.type = ele_type;
    tmp.ele = &ele_list[ele_type];
    tmp.pos = pos;
    tmp.force = pos*0;
    tmp.if_move = 2;
    atm_list.push_back(tmp);
    num_atm++;
    count_move_atoms();
}

void cell :: rm_atom(int ind_atm) {
    if (atm_list[ind_atm].if_move <= 1) {
        cout<<"Error: Can not remove atom "<<ind_atm+1<<' '<<atm_list[ind_atm].ele->sym<<", not removable"<<endl;
        exit(EXIT_FAILURE);
    }
    atm_list.erase(atm_list.begin() + ind_atm);
    num_atm--;
    count_move_atoms();
}

void cell :: sp_atom(int s1, int s2) {
    vec tmp;
    tmp = atm_list[s1].pos;
    atm_list[s1].pos = atm_list[s2].pos;
    atm_list[s2].pos = tmp;

    tmp = atm_list[s1].force;
    atm_list[s1].force = atm_list[s2].force;
    atm_list[s2].force = tmp;
}

void cell :: update_tb(double T) {
    for(int t1=0; t1<num_ele; t1++)
        ele_list[t1].update_tb(T);
}

void cell :: min_distance(vec pos, double& rr, int& ind) {
    double r_tmp;
    rr = 1e10;
    for(int t1=-1; t1<2; t1++)
        for(int t2=-1; t2<2; t2++)
            for(int t3=-1; t3<2; t3++)
                for(int t4=0; t4<num_atm; t4++) {
                    r_tmp=(pos + latt[0]*t1 + latt[1]*t2 + latt[2]*t3 - atm_list[t4].pos).norm();
                    if (rr > r_tmp) {
                        rr = r_tmp;
                        ind = t4;
                    }
                }
}

void cell :: print() {
    cout<<"Number of elements: "<<num_ele<<endl;
    cout<<"Number of atoms: "<<num_atm<<endl;
    cout<<"Threshold of height: "<<h_min<<", "<<h_max<<endl;
    cout<<"Energy: "<<energy<<endl;
    cout<<"Volume: "<<vol<<endl;
    cout<<"Total removable atoms: "<<num_atm_remove<<endl;
    cout<<"Total movable atoms: "<<num_atm_move<<endl;
    cout<<"Atoms per elements: "<<endl;
    for(int t1=0; t1<num_ele; t1++)
        cout<<ele_list[t1].sym<<'\t'<<num_ele_each[t1]<<endl;
    cout<<"Reovable atoms per elements: "<<endl;
    for(int t1=0; t1<num_ele; t1++)
        cout<<ele_list[t1].sym<<'\t'<<num_ele_each_remove[t1]<<endl;
    cout<<"Movable atoms per elements: "<<endl;
    for(int t1=0; t1<num_ele; t1++)
        cout<<ele_list[t1].sym<<'\t'<<num_ele_each_move[t1]<<endl;
//	cout<<"List of elements:"<<endl;
//	for(int t1=0; t1<num_ele; t1++)
//		ele_list[t1].print();
    cout<<"List of atoms:"<<endl;
    for(int t1=0; t1<num_atm; t1++)
        atm_list[t1].print();
}
