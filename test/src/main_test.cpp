//
// Created by leonard on 08.09.20.
//
#include <iostream>

#include "sample_test.h"
#include "common_test.h"

int main() {

    if (!TEST_ALL_COMMON()) {
        std::cerr << "Failed vector/poly test" << std::endl;
    }

    if (!TEST_ALL_SAMPLE()) {
        std::cerr << "Failed sample test" << std::endl;
    }

}
