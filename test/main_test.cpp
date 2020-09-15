//
// Created by leonard on 08.09.20.
//
#include <iostream>
#include "common.h"
#include "samples.h"

#include <chrono>

int main() {
    std::cout << "Hello World" << std::endl;

    NTT_engine<uint64_t, 1024, 134215681> engine{};
    TemplateVector<uint64_t, 1024, 134215681> vec(1u);

    TemplatePolynomial<uint64_t , 1024, 134215681> beef(engine, vec);
    TemplateRGswSample<uint32_t, 1024, 134215681> t1;

    auto t2 = t1 + t1;


    auto start = std::chrono::high_resolution_clock::now();
    auto beef1 = beef * beef;
    auto stop = std::chrono::high_resolution_clock::now();

    auto el = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

    std::cout << beef1[0] << " " << el.count();
}
