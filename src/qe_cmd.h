#ifndef QE_CMD
#define QE_CMD

#include <fstream>
#include <string>
#include "class.h"
#include "cell.h"

class qe_cmd
{
private:
	int num_core;
	int npool;
	int ndiag;
	std::string qe_exe;
	std::string mpi_launcher;

public:
	void read_from_in(std::ifstream& in);
	void write_qe_in(std::ifstream& in, std::ofstream& out, cell& c_new);
	void call(int if_test);
};

#endif
