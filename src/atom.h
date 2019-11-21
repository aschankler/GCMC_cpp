#ifndef ATOM
#define ATOM

#include <string>
#include <sstream>
#include <vector>
#include "class.h"
#include "element.h"
#include "vec.h"

class atom
{
private:
	// index of type to atom
	int type;
	// pointer to the element
	element *ele;
	// position
	vec pos;
	// force
	vec force;
	// define whether movable, 0, not movable; 1, movable but not removable; 2, all free
	int if_move;
friend cell;
friend mc;
friend qe_cmd;
public:
	void line_from_in(std::stringstream& ss, vector<element>& ele_list);
	void line_from_qe(std::stringstream& ss, vector<element>& ele_list);

	void print();
};
#endif
