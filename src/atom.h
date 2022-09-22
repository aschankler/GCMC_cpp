#ifndef __ATOM__
#define __ATOM__

#include <string>
#include <sstream>
#include <vector>
#include "element.h"
#include "vec.h"

class atom
{
public:
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

	void line_from_in(std::string tmp, std::vector<element>& ele_list);
	void line_from_qe(std::stringstream& ss, std::vector<element>& ele_list);

	void print() const;
};
#endif
