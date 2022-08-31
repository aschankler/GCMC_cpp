#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib>
#include <chrono>
#include "element.h"
#include "cell.h"
#include "mc.h"
#include "calculator.h"

using namespace std;

int main()
{
	chrono::high_resolution_clock::time_point now = chrono::high_resolution_clock::now();
	srand(now.time_since_epoch().count());
	// input and output files
	ifstream input;
	ofstream log, opt_axsf, accept_axsf, trial_axsf, trial_xsf, prop_axsf;
	// define of system
	cell sys_accept, sys_trial;
	mc mc_control;
	Calculator calculator_control;

	// read parameter file
	input.open("gcmc.in");
	sys_accept.read_from_in(input);
	mc_control.read_from_in(input);
	calculator_control.read_from_in(input);

	// initialize cell
	sys_accept.count_move_atoms();
	sys_accept.get_volume();
	sys_accept.update_tb(mc_control.temperature);

	// run first calculation
	cout<<"==============Begin iteration"<<setw(5)<<1<<"=============="<<endl;
	calculator_control.write_input(input,sys_accept);
	cout<<"Call calculator for the initial structure"<<endl;
	calculator_control.call(mc_control.if_test);
	if (!mc_control.if_test)
	{
		sys_accept.read_output(calculator_control.get_type());
		sys_accept.get_volume();
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
	log<<setw(9)<<"Iteration";
	for(int t1=0; t1<sys_accept.num_ele; t1++)
		log<<setw(4)<<sys_accept.ele_list[t1].sym;
	log<<setw(20)<<"DFT_E"<<setw(20)<<"Trial_E_f"<<setw(20)<<"Previous_E_f"<<setw(20)<<"Accept_E_f"<<setw(20)<<"Optimal_E_f"<<setw(14)<<"If_accept"<<endl;
	log<<setw(9)<<1;
	for(int t1=0; t1<sys_accept.num_ele; t1++)
		log<<setw(4)<<sys_accept.num_ele_each[t1];
	log<<fixed<<setprecision(9)<<setw(20)<<sys_accept.energy<<setw(20)<<mc_control.opt_e<<setw(20)<<mc_control.opt_e<<setw(20)<<mc_control.opt_e<<setw(20)<<mc_control.opt_e<<setw(14)<<1<<endl;
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
			sys_trial.get_volume();
		}
		else
		{
			sys_trial.energy += ((double)rand()/RAND_MAX - 0.5)*0.1;
		}
		// check if accept new structure
		if (mc_control.check_if_accept(sys_accept,sys_trial))
		{
			sys_accept = sys_trial;
			log<<setw(9)<<iter;
			for(int t1=0; t1<sys_trial.num_ele; t1++)
				log<<setw(4)<<sys_trial.num_ele_each[t1];
			log<<fixed<<setprecision(9)<<setw(20)<<sys_trial.energy<<setw(20)<<mc_control.e2<<setw(20)<<mc_control.e1<<setw(20)<<mc_control.e2<<setw(20)<<mc_control.opt_e<<setw(14)<<1<<endl;
		}
		else
		{
			log<<setw(9)<<iter;
			for(int t1=0; t1<sys_trial.num_ele; t1++)
				log<<setw(4)<<sys_trial.num_ele_each[t1];
			log<<fixed<<setprecision(9)<<setw(20)<<sys_trial.energy<<setw(20)<<mc_control.e2<<setw(20)<<mc_control.e1<<setw(20)<<mc_control.e1<<setw(20)<<mc_control.opt_e<<setw(14)<<0<<endl;
		}
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
