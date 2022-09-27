#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include "element.h"
#include "atom.h"
#include "vec.h"

using namespace std;

Atom atom_from_input(const std::string &tmp, vector<Element>& ele_list) {
	string ele_symbol;
	stringstream ss;

	Atom at;
	ss<<(tmp);
	ss>>ele_symbol>>at.pos_>>at.if_move_;
	at.type_ = -1;
	for (size_t t1=0 ; t1 < ele_list.size() ; t1++)
	{
		if (ele_symbol == ele_list[t1].sym_)
		{
			at.type_ = t1;
			break;
		}
	}
	if (at.type_ < 0)
	{
		cout<<"Error: Element "<<ele_symbol<<" not defined!"<<endl;
		exit(EXIT_FAILURE);
	}
	return at;
}

void Atom::print() const
{
	cout<<type_<<'\t'<<pos_<<'\t'<<if_move_<<endl;
}
