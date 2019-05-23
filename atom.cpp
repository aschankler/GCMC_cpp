#include <iostream>
#include <cstring>
#include <sstream>
#include "atom.h"
#include "vec.h"
#include "elements.h"

using namespace std;

atom :: atom()
{
	type = -1;
	pos.clean();
	force.clean();
	if_remove = 1;
	next = nullptr;
}

void atom :: get_type(string sym, elements& el)
{
	for(int t1=0; t1<el.num_el; t1++)
		if(sym == el.sym[t1])
		{
			type = t1;
			break;
		}
}

void atom :: line_from_in(string& info, elements& el)
{
	stringstream ss;
	string tmp;

	ss << (info);
	ss>>tmp>>pos>>if_remove;
	this->get_type(tmp,el);
}

void atom :: print()
{
	cout<<type<<'\t'<<pos<<'\t'<<if_remove<<endl;
}
