#include <iostream>
#include <fstream>
#include <cstring>
#include "element.h"
#include "cell.h"

using namespace std;

int main()
{
	ifstream input;
	input.open("in.dat");
	cell sys1;
	
	sys1.read_from_in(input);
	sys1.print();
	return 0;
}
