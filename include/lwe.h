//
// Created by leonard on 08.09.20.
//

#ifndef LWE_STRUCTS_LWE_H
#define LWE_STRUCTS_LWE_H

#include <cstddef>

#include "common.h"
#include "samples.h"
#include "random_base.h"

#define NOISE

enum KeyType {
    BINARY,
    ARBITRARY
};

/* for clarity */
template<typename T, size_t len, T modulus>
struct TemplateLweKey : TemplateVector<T, len, modulus> {
    KeyType type = ARBITRARY;
};

template<typename T, size_t dim, T modulus>
struct TemplateRingKey : TemplatePolynomial<T, dim, modulus> {
    KeyType type = ARBITRARY;
};

template<typename T, size_t vector_len, size_t poly_dim, T lwe_modulus, T ring_modulus>
struct TemplateCryptoEngine {

    /*
    TemplateCryptoEngine(double normal_mean, double normal_stddev, T base, T dim) : random_engine(lwe_modulus, ring_modulus, normal_mean, normal_stddev), ntt_engine(), rgsw_base(base), rgsw_decomp_dim(dim) {

    }
     */

    explicit TemplateCryptoEngine(RandomEngine<T>& base_engine, T base, T dim) : rgsw_base(base), rgsw_decomp_dim(dim), random_engine(base_engine), ntt_engine(){};

    TemplateCryptoEngine(const TemplateCryptoEngine& other) : random_engine(other.random_engine), ntt_engine(), rgsw_base(other.rgsw_base), rgsw_decomp_dim(other.rgsw_decomp_dim)  {

    }

    /*
     * @brief Encrypts a value \message by creating an LWE sample (a, <a,s> + e) and adding (0, \message)
     */
    TemplateLweSample<T, vector_len, lwe_modulus> encrypt_lwe(const TemplateLweKey<T, vector_len, lwe_modulus>& key, T message) {
        TemplateVector<T, vector_len, lwe_modulus> ct;

        for(int i = 0; i < vector_len; i++) {
            ct[i] = random_engine.generate_uniform_lwe();
        }

        auto b = ct * key;
        b = (b + message) % lwe_modulus;
#ifdef NOISE
        b = (b + random_engine.generate_normal()) % lwe_modulus;
#endif
        return TemplateLweSample<T, vector_len, lwe_modulus>(ct, b);
    }

    /*
     * @brief Encrypts a polynomial \message by creating a RLWE sample (a, a * s + e) and adding (0, \message)
     */
    TemplateRLweSample<T, poly_dim, ring_modulus> encrypt_rlwe(const TemplateRingKey<T, poly_dim, ring_modulus>& key, const TemplatePolynomial<T, poly_dim, ring_modulus>& message) {
        TemplatePolynomial<T, poly_dim, ring_modulus> ct_a(ntt_engine);

        for(int i = 0; i < poly_dim; i++) {
            ct_a[i] = random_engine.generate_uniform_rlwe();
        }

        auto ct_b = ct_a * key;
        ct_b += message;
#ifdef NOISE
        TemplatePolynomial<T, poly_dim, ring_modulus> noise(ntt_engine);
        for(int i = 0; i < poly_dim; i++) {
            noise[i] = random_engine.generate_normal();
        }
        ct_b += noise;
#endif
        return TemplateRLweSample<T, poly_dim, ring_modulus>(ct_a, ct_b);
    }

    /*
     * @brief Encrypts a polynomial \message by creating a RLWE' sample (vector of \rgsw_decomp_dim RLWE samples) and adding
     * a vector of \rgsw_decomp_dim (0, m_i) tuples, where m_i = \message * \rgsw_base ^ i
     */
    TemplateRLwePrimeSample<T, poly_dim, ring_modulus> encrypt_rlwe_prime(const TemplateRingKey<T, poly_dim, ring_modulus>& key, TemplatePolynomial<T, poly_dim, ring_modulus> message) {

        TemplateRLwePrimeSample<T, poly_dim, ring_modulus> ret(rgsw_base, rgsw_decomp_dim);
        for(int i = 0; i < rgsw_decomp_dim; i++) {
            ret[i] = encrypt_rlwe(key, message);
            message *= rgsw_base;
        }

        return ret;
    }

    /*
     * @brief Encrypts a polynomial \message by creating a RGSW sample where the first RLWE' sample encrypts (-\key * \message)
     * and the second one encrypts \message
     */
    TemplateRGswSample<T, poly_dim, ring_modulus> encrypt_rgsw(const TemplateRingKey<T, poly_dim, ring_modulus>& key, const TemplatePolynomial<T, poly_dim, ring_modulus>& message) {

        TemplateRGswSample<T, poly_dim, ring_modulus> ret(rgsw_base, rgsw_decomp_dim);

        auto sm = key * message;
        sm.negate();

        ret[0] = encrypt_rlwe_prime(key, sm);
        ret[1] = encrypt_rlwe_prime(key, message);

        return ret;
    }

    TemplateRGswSample<T, poly_dim, ring_modulus> trivial_rgsw(const TemplatePolynomial<T, poly_dim, ring_modulus>& msg) {
        return TemplateRGswSample<T, poly_dim, ring_modulus>(msg, rgsw_base, rgsw_decomp_dim);
    }

    T decrypt_lwe(const TemplateLweSample<T, vector_len, lwe_modulus>& ct, const TemplateLweSample<T, vector_len, lwe_modulus>& key) {
        return ct.getB() - ct.getA() * key;
    }

    TemplatePolynomial<T, poly_dim, ring_modulus> decrypt_rlwe(const TemplateRLweSample<T, poly_dim, ring_modulus>& ct, const TemplateRingKey<T, poly_dim, ring_modulus>& key) {
        auto pt = ct.getB() - ct.getA() * key;
        pt.setFormat(DEFAULT);
        return pt;
    }

    TemplatePolynomial<T, poly_dim, ring_modulus> decrypt_rlwe_prime(const TemplateRLwePrimeSample<T, poly_dim, ring_modulus>& ct, const TemplateRingKey<T, poly_dim, ring_modulus>& key) {
        return decrypt_rlwe(ct[0], key);
    }

    TemplatePolynomial<T, poly_dim, ring_modulus> decrypt_rgsw(const TemplateRGswSample<T, poly_dim, ring_modulus>& ct, const TemplateRingKey<T, poly_dim, ring_modulus>& key) {
        return decrypt_rlwe_prime(ct.getRHS(), key);
    }

    TemplateRingKey<T, poly_dim, ring_modulus> generate_ring_key(KeyType type = ARBITRARY) {
        TemplateRingKey<T, poly_dim, ring_modulus> ret;

        if(type == ARBITRARY) {
            for (int i = 0; i < poly_dim; i++) {
                ret[i] = random_engine.generate_uniform_rlwe();
            }
        } else {
            for(int i = 0; i < poly_dim; i++) {
                ret[i] = random_engine.generate_uniform_binary();
            }
        }

        return ret;
    }

    TemplateLweKey<T, vector_len, lwe_modulus> generate_lwe_key(KeyType type = ARBITRARY) {
        TemplateLweKey<T, vector_len, lwe_modulus> ret;

        if (type == ARBITRARY) {
            for (int i = 0; i < vector_len; i++) {
                ret[i] = random_engine.generate_uniform_lwe();
            }
        } else {
            for (int i = 0; i < vector_len; i++) {
                ret[i] = random_engine.generate_uniform_binary();
            }
        }

        return ret;
    }

private:
    T rgsw_base;
    T rgsw_decomp_dim;

    RandomEngine<T>& random_engine;
    NTT_engine<T, poly_dim, ring_modulus> ntt_engine;
};

#endif //LWE_STRUCTS_LWE_H
