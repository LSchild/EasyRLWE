//
// Created by leonard on 16.09.20.
//


#include "common.h"
#include "parameters.h"

bool TEST_VECTOR_OPS() {

    Vector a,b;
    ELEMENT_TYPE x = 0;

    for(int i = 0; i < LWE_DIM; i++) {
        a[i] = i;
        b[i] = LWE_MOD - i;
    }

    auto ab = a;

    /* sum test */
    auto c = a + b;
    a += b;

    for(int i = 0; i < LWE_DIM; i++) {
        if(c[i] != a[i]) {
            std::cerr << __FILE__ << " " << __LINE__ << " ] Failed (a+=b) == (c = a + b)" << std::endl;
            return false;
        } else if(c[i] != 0) {
            std::cerr << __FILE__ << " " << __LINE__ << "Failed modulo, expected 0 got " << c[i] << std::endl;
            return false;
        }
    }

    /* diff test */
    c = a - b;
    a -= b;

    for(int i = 0; i < LWE_DIM; i++) {
        if(c[i] != a[i]) {
            std::cerr << __FILE__ << " " << __LINE__ << " ] Failed (a-=b) == (c = a - b)" << std::endl;
            return false;
        } else if(c[i] != a[i]) {
            std::cerr << __FILE__ << " " << __LINE__ << " Failed modulo, expected 0 got " << c[i] << std::endl;
            return false;
        }
    }

    /* scale test */
    c = a * 2u;
    a *= 2u;

    for(int i = 0; i < LWE_DIM; i++) {
        if(c[i] != a[i]) {
            std::cerr << __FILE__ << " " << __LINE__ << " ] Failed (a*=2) == (c = a * 2)" << std::endl;
            return false;
        } else if(c[i] != (x = mod_mul<ELEMENT_TYPE, RLWE_MOD>(2, ab[i]))) {
            std::cerr << __FILE__ << " " << __LINE__ << " Failed modulo, expected " << x << " got " << c[i] << std::endl;
            return false;
        }
    }

    /* dot product test */
    auto dotp = ab * ab;
    if(dotp != 768) {
        std::cerr << __FILE__ << " " << __LINE__ << " ] Failed dot product test. Expected 768, got " << dotp << std::endl;
        return false;
    }

    return true;
}

bool TEST_POLY_OPS() {

}

bool TEST_ALL_COMMON() {
    bool ret = TEST_VECTOR_OPS();

    return ret;
}