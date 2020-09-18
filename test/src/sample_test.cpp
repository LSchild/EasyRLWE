//
// Created by leonard on 18.09.20.
//

#include "ntt.h"
#include "random_base.h"
#include "lwe.h"
#include "samples.h"
#include "parameters.h"

#include "sample_test.h"

bool TEST_RGSW_MUL(CryptoEngine& cryptoEngine, NTT_engine<ELEMENT_TYPE, RLWE_DIM, RLWE_MOD>& nttEngine) {

    Polynomial a(nttEngine);
    a[0] = 10;

    auto ring_key = cryptoEngine.generate_ring_key();
    auto rlwe_ct = cryptoEngine.encrypt_rlwe(ring_key, a);
    auto rgsw_ct = cryptoEngine.encrypt_rgsw(ring_key, a);

    auto product = rgsw_ct * rlwe_ct;
    auto message = cryptoEngine.decrypt_rlwe(product, ring_key);

    std::cerr << message[0] << std::endl;
    return true;
}

bool TEST_ALL_SAMPLE() {

    /* setup */
    NTT_engine<ELEMENT_TYPE, RLWE_DIM, RLWE_MOD> nttEngine;
    DefaultRandomEngine<ELEMENT_TYPE> randomEngine(LWE_MOD, RLWE_MOD, 0, 3.19);
    CryptoEngine cryptoEngine(randomEngine, 1024u,  3);

    TEST_RGSW_MUL(cryptoEngine, nttEngine);

    return true;
}