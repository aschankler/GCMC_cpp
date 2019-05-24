#include <iostream>
#include <fstream>
#include <cstring>
#include "element.h"
#include "cell.h"
#include "mc.h"

using namespace std;

int main()
{
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
	sys2 = sys1;
	sys2.read_from_qe(qe_out);
	sys2.ad_atom(pos1,0);
	sys2.ad_atom(pos1,1);
	sys2.ad_atom(pos1,2);
	sys2.ad_atom(pos1,2);
	sys2.rm_atom(4);
	sys2.sp_atom(0,6);
	sys2.print();
	cout<<endl;
	run1.print();

	for(size_t t1=0; t1<10; t1++)
	{
		sys2.min_distance(pos1*(t1/5.0),rr,ind);
	}
	return 0;
}
