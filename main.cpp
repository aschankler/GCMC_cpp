#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include "element.h"
#include "cell.h"
#include "mc.h"

using namespace std;

int main()
{
	srand(time(0));
	ifstream input, qe_out;
	input.open("in.dat");
	qe_out.open("t.out");
	cell sys1, sys2;
	mc run1;
	double pos0[3]={1.1,2.2,3.3}, rr;
	int ind;
	vec pos1;
	pos1=&pos0[0];
	
	run1.read_from_in(input);
	sys1.read_from_in(input);
	sys1.count_move_atoms();
	sys1.get_volume();
	sys1.update_tb(run1.T);
	sys1.read_from_qe(qe_out);
	sys1.ad_atom(pos1,0);

	sys2 = sys1;
	for (size_t t1=0; t1<10; t1++)
	{
		run1.create_new_structure(sys2,sys2);
		cout<<"----------------"<<endl;
	}
	/*
	sys2.rm_atom(4);
	sys2.sp_atom(0,6);
	*/
	sys1.print();
	cout<<"----------------"<<endl;
	sys2.print();
	return 0;
}
