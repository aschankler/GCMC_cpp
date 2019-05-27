#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include "element.h"
#include "cell.h"
#include "mc.h"
#include "qe_cmd.h"

using namespace std;

int main()
{
	srand(time(0));
	// input and output files
	ifstream input, qe_out;
	ofstream qe_in, log, opt_axsf, accept_axsf, trial_axsf;
	// define of system
	cell sys_accept, sys_trial;
	mc mc_control;
	qe_cmd qe_control;

	// read parameter file
	input.open("param.in");
	sys_accept.read_from_in(input);
	mc_control.read_from_in(input);
	qe_control.read_from_in(input);

	// initialize cell
	sys_accept.count_move_atoms();
	sys_accept.get_volume();
	sys_accept.update_tb(mc_control.T);

	// run first QE
	cout<<"==============Begin iteration"<<setw(5)<<1<<"=============="<<endl;
	qe_in.open("qe.in"); qe_control.write_qe_in(input,qe_in,sys_accept); qe_in.close();
	qe_control.call(mc_control.if_test);
	qe_out.open("qe.out"); sys_accept.read_from_qe(qe_out); qe_out.close();
	opt_axsf.open("save_opt.axsf");
	accept_axsf.open("save_accept.axsf");
	trial_axsf.open("save_trial.axsf");
	sys_accept.write_axsf(opt_axsf);
	sys_accept.write_axsf(accept_axsf);
	sys_accept.write_axsf(trial_axsf);
	sys_accept.write_axsf(opt_axsf,1);
	sys_accept.write_axsf(accept_axsf,1);
	sys_accept.write_axsf(trial_axsf,1);

	// initialize mc
	mc_control.save_opt_structure(sys_accept);

	// start mc iteration
	log.open("log.dat");
	log<<setw(9)<<"Iteration"<<setw(20)<<"DFT_E"<<setw(20)<<"Trial_E_f"<<setw(20)<<"Previous_E_f"<<setw(20)<<"Accept_E_f"<<setw(20)<<"Optimal_E_f"<<setw(14)<<"If_accept"<<endl;
	log<<setw(9)<<1<<fixed<<setprecision(9)<<setw(20)<<sys_accept.energy<<setw(20)<<mc_control.opt_e<<setw(20)<<mc_control.opt_e<<setw(20)<<mc_control.opt_e<<setw(20)<<mc_control.opt_e<<setw(14)<<1<<endl;
	cout<<"================================================"<<endl<<endl;
	for(size_t iter=2; iter<=mc_control.max_iter; iter++)
	{
		cout<<"==============Begin iteration"<<setw(5)<<iter<<"=============="<<endl;
		// create new structure
		mc_control.create_new_structure(sys_accept,sys_trial);
		qe_in.open("qe.in"); qe_control.write_qe_in(input,qe_in,sys_trial); qe_in.close();
		// run QE
		qe_control.call(mc_control.if_test);
		// read QE result
		if (!mc_control.if_test)
		{
			qe_out.open("qe.out"); sys_trial.read_from_qe(qe_out); qe_out.close();
		}
		else
		{
			sys_trial.energy += ((double)rand()/RAND_MAX - 0.5)*0.0;
		}
		// check if accept new structure
		if (mc_control.check_if_accept(sys_accept,sys_trial))
		{
			sys_accept = sys_trial;
			log<<setw(9)<<iter<<fixed<<setprecision(9)<<setw(20)<<sys_trial.energy<<setw(20)<<mc_control.e2<<setw(20)<<mc_control.e1<<setw(20)<<mc_control.e2<<setw(20)<<mc_control.opt_e<<setw(14)<<1<<endl;
		}
		else
		{
			log<<setw(9)<<iter<<fixed<<setprecision(9)<<setw(20)<<sys_trial.energy<<setw(20)<<mc_control.e2<<setw(20)<<mc_control.e1<<setw(20)<<mc_control.e1<<setw(20)<<mc_control.opt_e<<setw(14)<<0<<endl;
		}
		// write to axsf file
		mc_control.opt_c.write_axsf(opt_axsf,iter);
		sys_accept.write_axsf(accept_axsf,iter);
		sys_trial.write_axsf(trial_axsf,iter);
		cout<<"================================================"<<endl<<endl;
	}
	return 0;
}
