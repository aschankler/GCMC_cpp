#ifndef __ATOM__
#define __ATOM__

#include <string>
#include <vector>

#include "element.h"
#include "vec.h"

class Atom {
public:
	// index of type to atom
	int type;
	// pointer to the element
	Element *ele;
	// position
	vec pos;
	// force
	vec force;
	// define whether movable, 0, not movable; 1, movable but not removable; 2, all free
	int if_move;

	void line_from_in(std::string tmp, std::vector<Element>& ele_list);

	void print() const;
};
#endif
