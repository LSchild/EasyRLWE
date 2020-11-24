//
// Created by leonard on 08.09.20.
//
#include <iostream>

#include "sample_test.h"
#include "common_test.h"

int main() {

    bool all = true;
    bool ret = TEST_ALL_COMMON();

    all &= ret;
    if (!ret) {
        std::cerr << "Failed vector/poly test" << std::endl;
    }

    ret = TEST_ALL_SAMPLE();
    all &= ret;
    if (!ret) {
        std::cerr << "Failed sample test" << std::endl;
    }


    if(all) {
        std::cerr << "All tests passed" << std::endl;
    }

}
