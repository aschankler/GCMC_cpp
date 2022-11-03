#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <string>
#include <cstdlib>
#include "element.h"
#include "cell.h"
#include "mc.h"
#include "rng.h"
#include "calculator.h"

using namespace std;

namespace gcmc {
	std::string version = "0.3.0";
}


// From: https://stackoverflow.com/a/868894
class CliParser {
    public:
        CliParser(int &argc, char **argv) {
            for (int i = 1; i < argc; ++i)
                tokens_.push_back(std::string(argv[i]));
        }
        const std::string& get_option(const std::string &option) const {
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(tokens_.begin(), tokens_.end(), option);
            if (itr != tokens_.end() && ++itr != tokens_.end()) {
                return *itr;
            }
            static const std::string empty_string("");
            return empty_string;
        }
        bool option_exists(const std::string &option) const {
            return std::find(tokens_.begin(), tokens_.end(), option) != tokens_.end();
        }
    private:
        std::vector<std::string> tokens_;
};


int main(int argc, char *argv[])
{
	// Parse CLI options
	CliParser options(argc, argv);

	if (options.option_exists("-v") || options.option_exists("--version")) {
		cout << "aiGCMC " << gcmc::version << endl;
		exit(EXIT_SUCCESS);
	}

	std::string in_file;
	if (options.option_exists("-in")) {
		in_file = options.get_option("-in");
	} else {
		in_file = "gcmc.in";
	}

	// input and output files
	ofstream log, opt_axsf, accept_axsf, trial_axsf, trial_xsf, prop_axsf;

	gcmc::init_rng();

	// read parameter file
	ifstream input(in_file);
	Cell cell_accept = cell_from_in(input);
	MCMC mc_control = mcmc_from_in(input);
	Calculator calculator = calculator_from_in(input);

	// initialize cell
	cell_accept.update_tb(mc_control.temperature);

	// run first calculation
	cout<<"==============Begin iteration"<<setw(5)<<1<<"=============="<<endl;
	calculator.write_input(input, cell_accept);
	cout<<"Call calculator for the initial structure"<<endl;
	calculator.call(mc_control.if_test);

	if (!mc_control.if_test) {
		cell_accept.read_output(calculator.get_type());
	} else {
		cell_accept.energy = 0;
	}
	// print info
	cout<<endl<<"Save structures to files"<<endl;
	if (cell_accept.control.if_vc_relax_)
		cout<<"    Warning: Lattice constant is relaxed, .axsf is meanless"<<endl;

	opt_axsf.open("save_opt.axsf");
	accept_axsf.open("save_accept.axsf");
	trial_axsf.open("save_trial.axsf");
	trial_xsf.open("save_trial.xsf");
	prop_axsf.open("save_propose.axsf");

	cell_accept.write_axsf(opt_axsf);
	cell_accept.write_axsf(accept_axsf);
	cell_accept.write_axsf(trial_axsf);
	cell_accept.write_axsf(prop_axsf);
	cell_accept.write_axsf(opt_axsf,1);
	cell_accept.write_axsf(accept_axsf,1);
	cell_accept.write_axsf(trial_axsf,1);
	cell_accept.write_xsf(trial_xsf,1);
	cell_accept.write_axsf(prop_axsf, 1);

	// initialize mc
	mc_control.save_opt_structure(cell_accept);

	// start mc iteration
	log.open("log.dat");
    mc_control.log_header(log, cell_accept);
    mc_control.log_step(log, cell_accept, 1);

	cout<<"================================================"<<endl<<endl;
	for(int iter=2; iter<=mc_control.max_iter; iter++)
	{
		cout<<"==============Begin iteration"<<setw(5)<<iter<<"=============="<<endl;
		// create new structure
		Cell cell_trial = mc_control.create_new_structure(cell_accept);
		// Save initial proposal
		cell_trial.write_axsf(prop_axsf, iter);
		// execute calculator
		calculator.write_input(input, cell_trial);
		calculator.call(mc_control.if_test);

		// read calculator result
		if (!mc_control.if_test) {
			cell_trial.read_output(calculator.get_type());
		} else {
			cell_trial.energy += (gcmc::rand_uniform() - 0.5)*0.1;
		}

		// check if accept new structure
		if (mc_control.check_if_accept(cell_accept, cell_trial)) {
			cell_accept = cell_trial;
		}
        mc_control.log_step(log, cell_trial, iter);

		// write to axsf file
		cout<<endl<<"Save structures to files"<<endl;
		if (cell_accept.control.if_vc_relax_)
			cout<<"    Warning: Lattice constant is relaxed, .axsf is meanless"<<endl;
		mc_control.opt_c.write_axsf(opt_axsf,iter);
		cell_accept.write_axsf(accept_axsf,iter);
		cell_trial.write_axsf(trial_axsf,iter);
		cell_trial.write_xsf(trial_xsf,iter);
		cout<<"================================================"<<endl<<endl;
	}
	return 0;
}
