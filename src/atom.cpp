#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include "element.h"
#include "atom.h"
#include "vec.h"

using namespace std;

void atom :: line_from_in(string tmp, vector<element>& ele_list)
{
	string ele_symbol;
	stringstream ss;

	ss<<(tmp);
	ss>>ele_symbol>>pos>>if_move;
	ele = nullptr;
	for (size_t t1=0 ; t1 < ele_list.size() ; t1++)
	{
		if (ele_symbol == ele_list[t1].sym)
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

void atom :: print()
{
	cout<<ele->sym<<'\t'<<pos<<'\t'<<if_move<<endl;
//	ele->print();
}
