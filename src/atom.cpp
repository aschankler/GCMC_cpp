#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include "element.h"
#include "atom.h"
#include "vec.h"

using namespace std;

void Atom::line_from_in(string tmp, vector<Element>& ele_list)
{
	string ele_symbol;
	stringstream ss;

	ss<<(tmp);
	ss>>ele_symbol>>pos>>if_move;
	ele = nullptr;
	for (size_t t1=0 ; t1 < ele_list.size() ; t1++)
	{
		if (ele_symbol == ele_list[t1].sym_)
		{
			type = t1;
			ele = &ele_list[t1];
			break;
		}
	}
	if (ele == nullptr)
	{
		cout<<"Error: Element "<<ele_symbol<<" not defined!"<<endl;
		exit(EXIT_FAILURE);
	}
}

void Atom::print() const
{
	cout<<ele->sym_<<'\t'<<pos<<'\t'<<if_move<<endl;
}
