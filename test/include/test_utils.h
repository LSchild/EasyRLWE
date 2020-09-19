//
// Created by leonard on 19.09.20.
//

#ifndef LWE_STRUCTS_TEST_UTILS_H
#define LWE_STRUCTS_TEST_UTILS_H

#include <iostream>

inline void error(const char* msg, const char* file, const int line) {
    std::cerr << file << " " << line << " ] " << msg << std::endl;
}

#endif //LWE_STRUCTS_TEST_UTILS_H
