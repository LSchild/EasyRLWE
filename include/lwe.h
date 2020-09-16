//
// Created by leonard on 08.09.20.
//

#ifndef LWE_STRUCTS_LWE_H
#define LWE_STRUCTS_LWE_H

#include <cstddef>

#include "common.h"
#include "samples.h"
#include "random_base.h"

/* for clarity */
template<typename T, size_t len, T modulus>
using TemplateLweKey = TemplateVector<T, len, modulus>;

template<typename T, size_t dim, T modulus>
using TemplateRingKey = TemplatePolynomial<T, dim, modulus>;

template<typename T, size_t vector_len, size_t poly_dim, T lwe_modulus, T ring_modulus>
struct TemplateCryptoEngine {

    TemplateCryptoEngine(double normal_mean, double normal_stddev) : random_engine(lwe_modulus, ring_modulus, normal_mean, normal_stddev), ntt_engine() {
    }

    explicit TemplateCryptoEngine(const RandomEngine<T>& base_engine) : random_engine(base_engine), ntt_engine() {};

    TemplateCryptoEngine(const TemplateCryptoEngine& other) : random_engine(other.random_engine), ntt_engine()  {

    }

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

    TemplateRLweSample<T, poly_dim, ring_modulus> encrypt_rlwe(const TemplateRingKey<T, poly_dim, ring_modulus>& key, const TemplatePolynomial<T, poly_dim, ring_modulus>& message) const {
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
        return TemplatePolynomial<T, poly_dim, ring_modulus>(ct_a, ct_b);
    }

    TemplateRLwePrimeSample<T, poly_dim, ring_modulus> encrypt_rlwe_prime(const TemplateRingKey<T, poly_dim, ring_modulus>& key, TemplatePolynomial<T, poly_dim, ring_modulus> message) {

    }

private:
    RandomEngine<T> random_engine;
    NTT_engine<T, poly_dim, ring_modulus> ntt_engine;
};

#endif //LWE_STRUCTS_LWE_H
