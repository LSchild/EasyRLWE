//
// Created by leonard on 16.09.20.
//


#include "common.h"
#include "parameters.h"

#include <chrono>

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

    Polynomial a;
    NTT_engine<ELEMENT_TYPE, RLWE_DIM, RLWE_MOD> ntt;

    for(int i = 0; i < RLWE_DIM; i++) {
        a[i] = i;
    }

    /* constructors */
    Polynomial d(ntt);
    Polynomial e(a);
    Polynomial f(ntt, 1, 0);

    for(int i = 0; i < RLWE_DIM; i++) {
        d[i] = i;
    }

    /* equality */
    if(a != e) {
        std::cerr << __FILE__ << " " << __LINE__ << " ] Polynomials should match..." << std::endl;
        return false;
    }

    /* ntt */
    e.setFormat(NTT);
    e.setFormat(DEFAULT);

    if(a != e) {
        std::cerr << __FILE__ << " " << __LINE__ << " ] 2 * NTT should be identical " << std::endl;
        return false;
    }

    auto start = std::chrono::high_resolution_clock::now();

    f = d + d;

    auto stop = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    std::cerr << "Elapsed " << elapsed.count() << std::endl;

    for(int i = 0; i < RLWE_DIM; i++) {
        if(f[i] != 2 * i) {
            std::cerr << __FILE__ << " " << __LINE__ << " ] Sum does not match " << std::endl;
            return false;
        }
    }

    f -= d;
    for(int i = 0; i < RLWE_DIM; i++) {
        if(f[i] != d[i]) {
            std::cerr << __FILE__ << " " << __LINE__ << " ] 2 * d - d != d " << std::endl;
            return false;
        }
    }

    return true;
}

bool TEST_ALL_COMMON() {

    bool ret = TEST_VECTOR_OPS();
    ret &= TEST_POLY_OPS();

    return ret;
}