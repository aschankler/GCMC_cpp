#ifndef ATOM
#define ATOM

#include <cstring>
#include "class.h"
#include "vec.h"

class atom
{
private:
	int type;
	vec pos;
	vec force;
	int if_remove;
	atom *next;
friend cell;
public:
	atom();
	void get_type(std::string sym, elements& el);
	void line_from_in(std::string& info, elements& el);
	void line_from_qe(std::string& info, elements& el);

	void print();
};
#endif
