#ifndef __CALCULATOR__
#define __CALCULATOR__

#include <fstream>
#include <string>
#include "cell.h"

class calculator
{
public:
	int num_core;
	int npool;
	int ndiag;
	int calculator_type; // 1: qe, 2: vasp
	std::string exe;
	std::string mpi_launcher;

	void read_from_in(std::ifstream& in);
	void write_input(std::ifstream& in, cell& c_new);
	void write_qe_in(std::ifstream& in, cell& c_new);
	void write_vasp_in(std::ifstream& in, cell& c_new);
	void call(int if_test);
};

#endif
