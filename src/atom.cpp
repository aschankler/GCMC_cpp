#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include "atom.h"
#include "vec.h"
#include "element.h"

using namespace std;

void atom :: line_from_in(stringstream& ss, vector<element>& ele_list)
{
	string ele_symbol;
	ss>>ele_symbol>>pos>>if_move;
	auto num_ele = ele_list.size();
	for (size_t t1=0 ; t1 < num_ele ; t1++)
	{
		if (ele_symbol == ele_list[t1].sym)
		{
			type = t1;
			ele = &ele_list[t1];
			break;
		}
	}
}

void atom :: print()
{
	cout<<ele->sym<<'\t'<<pos<<'\t'<<if_move<<endl;
//	ele->print();
}
