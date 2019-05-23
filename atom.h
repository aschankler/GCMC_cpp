#ifndef ATOM
#define ATOM

#include <cstring>
#include <sstream>
#include <vector>
#include "class.h"
#include "element.h"
#include "vec.h"

class atom
{
private:
	int type;
	element *ele;
	vec pos;
	vec force;
	int if_remove;
friend cell;
public:
	void line_from_in(std::stringstream& ss, vector<element>& ele_list);
	void line_from_qe(std::stringstream& ss, vector<element>& ele_list);

	void print();
};
#endif
