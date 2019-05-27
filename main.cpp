#include <iostream>
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
	ifstream input, qe_out;
	ofstream qe_in, trial_axsf;
	cell sys1, sys2;
	mc run1;
	qe_cmd qe_control;
	double pos0[3]={1.1,2.2,3.3}, rr;
	int ind;
	vec pos1;
	pos1=&pos0[0];
	
	input.open("param.in");
	qe_in.open("qe.in");
	qe_out.open("qe.out");
	trial_axsf.open("trial.axsf");

	run1.read_from_in(input);
	sys1.read_from_in(input);
	qe_control.read_from_in(input);
	sys1.count_move_atoms();
	sys1.get_volume();
	sys1.update_tb(run1.T);
	sys1.read_from_qe(qe_out);
	sys1.ad_atom(pos1,0);

	run1.save_opt_structure(sys1);
	sys1.write_axsf(trial_axsf);
	for (size_t t1=0; t1<run1.max_iter; t1++)
	{
		run1.create_new_structure(sys1,sys2);
		sys2.energy += ((double)rand()/RAND_MAX - 0.5)*0.1;
		if (run1.if_accept(sys1,sys2))
			sys1 = sys2;
		cout<<"----------------"<<endl;
		sys2.write_axsf(trial_axsf,t1+1);
	}
	cout<<"----------------"<<endl;
	sys2.print();
	qe_control.write_qe_in(input,qe_in,sys2);
	//qe_control.call(run1.if_test);

	input.close();
	qe_in.close();
	qe_out.close();
	trial_axsf.close();
	return 0;
}
