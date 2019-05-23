#include <iostream>
#include <fstream>
#include <cstring>
#include "elements.h"
#include "cell.h"

using namespace std;

int main()
{
	ifstream input;
	input.open("in.dat");
	elements el_list;
	cell sys1;
	
	el_list.get_param(input);
	sys1.read_from_in(input, el_list);
	sys1.print();
	return 0;
}
