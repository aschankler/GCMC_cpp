#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "calculator.h"
#include "auxiliary.h"

using namespace std;


Calculator calculator_from_in(istream& in) {
    std::string exe_name;
    read(in, "exe", '=', exe_name);

    // Read parallelism
    int num_core;
    std::string mpi_launcher;
    read(in, "mpi_launcher", '=', mpi_launcher);
    read(in, "num_core", '=', num_core);

    std::string para_params;

    if (not read_opt(in, "para_params", '=', para_params)) {
        para_params = "";
    }

    // Invalid for vasp
    int npool;
    if (read_opt(in,"npool",'=',npool)) {
        para_params.append(" -nk " + std::to_string(npool));
    }
    int ndiag;
    if (read_opt(in,"ndiag",'=',ndiag)) {
        para_params.append(" -ndiag " + std::to_string(ndiag));
    }

    return Calculator(exe_name, mpi_launcher, num_core, para_params);
}


Calculator::Calculator(std::string name, std::string launcher, int nc, std::string params)
    : exe_name_(name), mpi_launcher_(launcher), num_core_(nc), para_params_(params) {
    if (exe_name_.find("pw") != string::npos) {
        calculator_type_ = 1;
    } else if(exe_name_.find("vasp") != string::npos) {
        calculator_type_ = 2;
    } else {
        throw std::runtime_error("Currently only support QE and VASP!");
    }
}

void Calculator::write_input(ifstream& in, const Cell& c_new) const {
    switch (calculator_type_) {
    case 1: { //QE
        write_qe_in(in,c_new);
        break;
    }
    case 2: { //vasp
        write_vasp_in(in,c_new);
        break;
    }
    default: {
        cout<<"Error: Invalid calculator type!"<<endl;
        exit(EXIT_FAILURE);
    }
    }
}

void Calculator::write_qe_in(ifstream& in, const Cell& c_new) const {
    ofstream out("qe.in");
    string label_qe_input = "begin_qe_input";
    string label_qe_input2 = "end_qe_input";
    string label_nat = "nat";
    string cell_param = "CELL_PARAMETERS (angstrom)";
    string position = "ATOMIC_POSITIONS (angstrom)";
    string tmp;

    in.clear();
    in.seekg(ios::beg);
    while(getline(in,tmp))
        if(tmp.find(label_qe_input) != string::npos)
            break;
    if(in.eof()) {
        cout<<"Error: Cannot find "<<label_qe_input<<'!'<<endl;
        exit(EXIT_FAILURE);
    }
    // copy parameter to qe.in
    while(getline(in,tmp)) {
        if(tmp.find(label_qe_input2) == string::npos) {
            if(tmp.find(label_nat) == string::npos)
                out<<tmp<<endl;
            else
                out<<tmp<<" "<<c_new.num_atm<<endl;
        } else
            break;
    }
    if(in.eof()) {
        cout<<"Error: Cannot find "<<label_qe_input2<<'!'<<endl;
        exit(EXIT_FAILURE);
    }
    // copy cell parameter and position to qe.in
    out<<endl<<cell_param<<endl;
    out<<c_new.latt[0]<<endl;
    out<<c_new.latt[1]<<endl;
    out<<c_new.latt[2]<<endl;

    out<<endl<<position<<endl;
    for(int t1=0; t1<c_new.num_atm; t1++) {
        out<<setw(2)<<c_new.atom_type(t1).sym_<<"    ";
        out<<c_new.atm_list[t1].pos_;
        if(c_new.atm_list[t1].if_move_ == 0)
            out<<" 0 0 0"<<endl;
        else
            out<<endl;
    }
    out.close();
}

void Calculator::write_vasp_in(std::ifstream &in, const Cell &c_new) const {
    ofstream incar("INCAR");
    ofstream kpoints("KPOINTS");
    ofstream poscar("POSCAR");
    string label_vasp_incar = "begin_vasp_incar";
    string label_vasp_incar2= "end_vasp_incar";
    string label_vasp_kpoints = "begin_vasp_kpoints";
    string label_vasp_kpoints2= "end_vasp_kpoints";
    string tmp;
    // generate incar
    in.clear();
    in.seekg(ios::beg);
    while(getline(in,tmp))
        if(tmp.find(label_vasp_incar) != string::npos)
            break;
    if(in.eof()) {
        cout<<"Error: Cannot find "<<label_vasp_incar<<'!'<<endl;
        exit(EXIT_FAILURE);
    }
    while(getline(in,tmp)) {
        if(tmp.find(label_vasp_incar2) == string::npos)
            incar<<tmp<<endl;
        else
            break;
    }
    if(in.eof()) {
        cout<<"Error: Cannot find "<<label_vasp_incar2<<'!'<<endl;
        exit(EXIT_FAILURE);
    }
    incar.close();
    // generate kpoint
    in.clear();
    in.seekg(ios::beg);
    while(getline(in,tmp))
        if(tmp.find(label_vasp_kpoints) != string::npos)
            break;
    if(in.eof()) {
        cout<<"Error: Cannot find "<<label_vasp_kpoints<<'!'<<endl;
        exit(EXIT_FAILURE);
    }
    while(getline(in,tmp)) {
        if(tmp.find(label_vasp_kpoints2) == string::npos)
            kpoints<<tmp<<endl;
        else
            break;
    }
    if(in.eof()) {
        cout<<"Error: Cannot find "<<label_vasp_kpoints2<<'!'<<endl;
        exit(EXIT_FAILURE);
    }
    kpoints.close();

    // generate poscar
    poscar<<"gcmc_auto_generated"<<endl;
    poscar<<"   1.000000000000"<<endl;
    poscar<<c_new.latt[0]<<endl;
    poscar<<c_new.latt[1]<<endl;
    poscar<<c_new.latt[2]<<endl;
    for(auto m1:c_new.ele_list)
        poscar<<setw(5)<<m1.sym_;
    poscar<<endl;
    for(auto m1:c_new.num_ele_each)
        poscar<<setw(5)<<m1;
    poscar<<endl;
    poscar<<"Selective dynamics"<<endl;
    poscar<<"Direct"<<endl;
    // atomic positions
    vec frac_coord;
    for(auto m1:c_new.ele_list) {
        for(auto m2:c_new.atm_list) {
            if(m1.sym_ == c_new.ele_list[m2.type_].sym_) {
                frac_coord.x[0] = m2.pos_*c_new.latt_inv[0];
                frac_coord.x[1] = m2.pos_*c_new.latt_inv[1];
                frac_coord.x[2] = m2.pos_*c_new.latt_inv[2];
                poscar<<frac_coord;
                if(m2.if_move_)
                    poscar<<"  T  T  T"<<endl;
                else
                    poscar<<"  F  F  F"<<endl;
            }
        }
    }
    poscar.close();
}


std::string Calculator::cli_string() const {
    stringstream ss;

    ss << mpi_launcher_ << " -n " << num_core_ << " " << exe_name_;

    switch(calculator_type_) {
    case 1: { // QE
        ss << " " << para_params_ << " -in qe.in > qe.out";
        break;
    }
    case 2: { // VASP
        ss << " > vasp.out";
        break;
    }
    default: {
        cout<<"Error: Invalid calculator type!"<<endl;
        exit(EXIT_FAILURE);
    }
    }

    return ss.str();
}


void Calculator::call(int if_test) const {
    string cli_str;

    cli_str = cli_string();
    cout << endl << "Launching calculator using command:" << endl;
    cout << "    " << cli_str << endl;

    if (not if_test) {
        system(cli_str.c_str());
        cout << "Calculation finished" << endl;
    }
}
