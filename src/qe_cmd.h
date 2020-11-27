#ifndef __QE_CMD__
#define __QE_CMD__

#include <fstream>
#include <string>
#include "cell.h"

class qe_cmd
{
public:
	int num_core;
	int npool;
	int ndiag;
	std::string qe_exe;
	std::string mpi_launcher;

	void read_from_in(std::ifstream& in);
	void write_qe_in(std::ifstream& in, std::ofstream& out, cell& c_new);
	void call(int if_test);
};

#endif
