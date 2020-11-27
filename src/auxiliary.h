#ifndef __AUXILIARY__
#define __AUXILIARY__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>

template <class T>
bool read(std::ifstream& in, std::string keyword, char dlmt, T& save)
{
	std::string tmp;
	std::stringstream ss;

	in.seekg(std::ios::beg); in.clear();
	while(getline(in,tmp))
	{
		if(tmp.find(keyword) != std::string::npos)
		{
			ss << (tmp);
			getline(ss,tmp,dlmt);
			ss >> save;
			return true;
		}
	}
	std::cout<<"Error: Can not find "<<keyword<<"!"<<std::endl;
	std::exit(EXIT_FAILURE);
	return false;
}
template <> inline
bool read<std::string>(std::ifstream& in, std::string keyword, char dlmt, std::string& save)
{
	std::string tmp;
	std::stringstream ss;

	in.seekg(std::ios::beg); in.clear();
	while(getline(in,tmp))
	{
		if(tmp.find(keyword) != std::string::npos)
		{
			ss << (tmp);
			getline(ss,tmp,dlmt);
			getline(ss,tmp,'"');
			getline(ss,save,'"');
			return true;
		}
	}
	std::cout<<"Error: Can not find "<<keyword<<"!"<<std::endl;
	std::exit(EXIT_FAILURE);
	return false;
}

#endif
