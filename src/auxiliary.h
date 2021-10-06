#ifndef __AUXILIARY__
#define __AUXILIARY__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>

/*
 * Search for `<keyword> <dlmt> <value>` lines and write the value to save
 * Return value indicates whether a matching line was found
 */
template <class T> bool read_opt(
    std::ifstream& in, std::string keyword, char dlmt, T& save
) {
    std::string line;
    std::stringstream ss;

    in.clear();
    in.seekg(0, in.beg);
	while(getline(in, line)) {
        if (line.find(keyword) != std::string::npos) {
            ss << (line);
            getline(ss, line, dlmt);
            ss >> save;
            return true;
        }
    }
    // No match found
    return false;
}

template <> inline bool read_opt<std::string>(
    std::ifstream& in, std::string keyword, char dlmt, std::string& save
) {
    std::string line;
    std::stringstream ss;

    in.clear();
    in.seekg(std::ios::beg);
    while(getline(in, line)) {
        if(line.find(keyword) != std::string::npos) {
            ss << (line);
            getline(ss,line,dlmt);
            getline(ss,line,'"');
            getline(ss,save,'"');
            return true;
        }
    }
    return false;
}

/*
 * Like `read_opt` but raises error if no match found
 */
template <class T> bool read(
    std::ifstream& in, std::string keyword, char dlmt, T& save
) {
    bool found = read_opt<T>(in, keyword, dlmt, save);
    if (not found) {
        std::cout<<"Error: Can not find "<<keyword<<"!"<<std::endl;
        std::exit(EXIT_FAILURE);
        return false;
    }
    return true;
}

#endif
