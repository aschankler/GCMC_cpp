#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib>
#include "element.h"
#include "cell.h"
#include "mc.h"
#include "rng.h"
#include "calculator.h"

using namespace std;


int main()
{
	// input and output files
	ifstream input;
	ofstream log, opt_axsf, accept_axsf, trial_axsf, trial_xsf, prop_axsf;
	// define of system
	Cell sys_accept, sys_trial;
	mc mc_control;
	Calculator calculator_control;

	gcmc::init_rng();

	// read parameter file
	input.open("gcmc.in");
	sys_accept.read_from_in(input);
	mc_control.read_from_in(input);
	calculator_control.read_from_in(input);

	// initialize cell
	sys_accept.update_tb(mc_control.temperature);

	// run first calculation
	cout<<"==============Begin iteration"<<setw(5)<<1<<"=============="<<endl;
	calculator_control.write_input(input,sys_accept);
	cout<<"Call calculator for the initial structure"<<endl;
	calculator_control.call(mc_control.if_test);
	if (!mc_control.if_test)
	{
		sys_accept.read_output(calculator_control.get_type());
	}
	else
	{
		sys_accept.energy = 0;
	}
	// print info
	cout<<endl<<"Save structures to files"<<endl;
	if (sys_accept.if_vc_relax)
		cout<<"    Warning: Lattice constant is relaxed, .axsf is meanless"<<endl;

	opt_axsf.open("save_opt.axsf");
	accept_axsf.open("save_accept.axsf");
	trial_axsf.open("save_trial.axsf");
	trial_xsf.open("save_trial.xsf");
	prop_axsf.open("save_propose.axsf");

	sys_accept.write_axsf(opt_axsf);
	sys_accept.write_axsf(accept_axsf);
	sys_accept.write_axsf(trial_axsf);
	sys_accept.write_axsf(prop_axsf);
	sys_accept.write_axsf(opt_axsf,1);
	sys_accept.write_axsf(accept_axsf,1);
	sys_accept.write_axsf(trial_axsf,1);
	sys_accept.write_xsf(trial_xsf,1);
	sys_accept.write_axsf(prop_axsf, 1);

	// initialize mc
	mc_control.save_opt_structure(sys_accept);

	// start mc iteration
	log.open("log.dat");
    mc_control.log_header(log, sys_accept);
    mc_control.log_step(log, sys_accept, 1);

	cout<<"================================================"<<endl<<endl;
	for(int iter=2; iter<=mc_control.max_iter; iter++)
	{
		cout<<"==============Begin iteration"<<setw(5)<<iter<<"=============="<<endl;
		// create new structure
		mc_control.create_new_structure(sys_accept,sys_trial);
		// Save initial proposal
		sys_trial.write_axsf(prop_axsf, iter);
		// execute calculator
		calculator_control.write_input(input,sys_trial);
		calculator_control.call(mc_control.if_test);
		// read calculator result
		if (!mc_control.if_test)
		{
			sys_trial.read_output(calculator_control.get_type());
		}
		else
		{
			sys_trial.energy += (gcmc::rand_uniform() - 0.5)*0.1;
		}

		// check if accept new structure
		if (mc_control.check_if_accept(sys_accept, sys_trial)) {
			sys_accept = sys_trial;
		}
        mc_control.log_step(log, sys_trial, iter);

		// write to axsf file
		cout<<endl<<"Save structures to files"<<endl;
		if (sys_accept.if_vc_relax)
			cout<<"    Warning: Lattice constant is relaxed, .axsf is meanless"<<endl;
		mc_control.opt_c.write_axsf(opt_axsf,iter);
		sys_accept.write_axsf(accept_axsf,iter);
		sys_trial.write_axsf(trial_axsf,iter);
		sys_trial.write_xsf(trial_xsf,iter);
		cout<<"================================================"<<endl<<endl;
	}
	return 0;
}
