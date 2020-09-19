//
// Created by leonard on 18.09.20.
//

#include "ntt.h"
#include "random_base.h"
#include "lwe.h"
#include "samples.h"
#include "parameters.h"
#include "test_utils.h"

#include "sample_test.h"

bool TEST_ENCRYPTION(CryptoEngine& cryptoEngine, NTT_engine<ELEMENT_TYPE, RLWE_DIM, RLWE_MOD>& nttEngine) {

    Polynomial message(nttEngine);
    for(int i = 0; i < RLWE_DIM; i++) {
        message[i] = i;
    }

    auto ring_key = cryptoEngine.generate_ring_key();

    auto rlwe = cryptoEngine.encrypt_rlwe(ring_key, message);
    auto rlwep = cryptoEngine.encrypt_rlwe_prime(ring_key, message);
    auto rgsw = cryptoEngine.encrypt_rgsw(ring_key, message);

    /* RLWE */
    auto dec = cryptoEngine.decrypt_rlwe(rlwe, ring_key);

    if(dec != message) {
        error("encryption and decryption do not match for RLWE", __FILE__, __LINE__);
        return false;
    }

    /* RLWE' */
    if(rlwep.getEntries().size() != rlwep.getDimension()) {
        error("dimension does not match actual length of entry vector", __FILE__, __LINE__);
        return false;
    }

    dec = cryptoEngine.decrypt_rlwe_prime(rlwep, ring_key);

    if(dec != message) {
        error("encryption and decryption do not match for RLWE", __FILE__, __LINE__);
        return false;
    }

    auto dec_1 = cryptoEngine.decrypt_rlwe(rlwep[1], ring_key);
    if(dec_1 != (message * rlwep.getBase())) {
        error("second slot of RLWE' should be equal to first times base", __FILE__, __LINE__);
        return false;
    }

    /* rgsw */
    dec = cryptoEngine.decrypt_rgsw(rgsw, ring_key);

    if(dec != message) {
        error("encryption and decryption do not match for RGSW", __FILE__, __LINE__);
        return false;
    }

    auto lhs = rgsw[0];
    auto dec_2 = cryptoEngine.decrypt_rlwe(lhs[0], ring_key);
    auto sm = ring_key * message;
    auto zero = dec_2 + sm;

    for(int i = 0; i < RLWE_DIM; i++) {
        if(zero[i] != 0) {
            error("s * m + lhs = s * m + (- s * m) should be equal to zero", __FILE__, __LINE__);
            return false;
        }
    }

    return true;
}

bool TEST_RGSW_MUL(CryptoEngine& cryptoEngine, NTT_engine<ELEMENT_TYPE, RLWE_DIM, RLWE_MOD>& nttEngine) {

    Polynomial message(nttEngine);
    message[0] = 32;

    auto ring_key = cryptoEngine.generate_ring_key();
    auto rgsw_ct = cryptoEngine.encrypt_rgsw(ring_key, message);
    auto rlwe_ct = cryptoEngine.encrypt_rlwe(ring_key, message);

    auto s = rgsw_ct * rlwe_ct;
    auto dec = cryptoEngine.decrypt_rlwe(s, ring_key);

    if(dec[0] != 1024) {
        error("RLWE times RGSW product does not work !", __FILE__, __LINE__);
        return false;
    }

    return true;
}

bool TEST_ALL_SAMPLE() {

    /* setup */
    NTT_engine<ELEMENT_TYPE, RLWE_DIM, RLWE_MOD> nttEngine;
    DefaultRandomEngine<ELEMENT_TYPE> randomEngine(LWE_MOD, RLWE_MOD, 0, 3.19);
    CryptoEngine cryptoEngine(randomEngine, 1024u,  3);

    TEST_RGSW_MUL(cryptoEngine, nttEngine);
    bool ret = TEST_ENCRYPTION(cryptoEngine, nttEngine);


    return ret;
}