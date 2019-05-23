#include <iostream>
#include <fstream>
#include <cstring>
#include "element.h"
#include "cell.h"

using namespace std;

int main()
{
	ifstream input, qe_out;
	input.open("in.dat");
	qe_out.open("t.out");
	cell sys1, sys2;
	
	sys1.read_from_in(input);
	sys2 = sys1;
	sys2.read_from_qe(qe_out);
	sys1.print();
	cout<<endl;
	sys2.print();
	return 0;
}
