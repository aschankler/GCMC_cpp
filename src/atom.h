#ifndef __ATOM__
#define __ATOM__

#include <string>
#include <vector>

#include "element.h"
#include "vec.h"

class Atom {
  public:
    // index of type to atom
    int type_;
    // position
    vec pos_;
    // force
    vec force_;
    // define whether movable, 0, not movable; 1, movable but not removable; 2, all free
    int if_move_;

    void print() const;
};

Atom atom_from_input(const std::string&, std::vector<Element>&);

#endif
