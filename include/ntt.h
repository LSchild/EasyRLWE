#ifndef NTT_H
#define NTT_H

#include <array>
#include <iostream>
#include <sys/mman.h>
#include <unistd.h>
#include <stdlib.h>
#include <immintrin.h>

#include "utils.h"

template<typename T, T modulus>
static T compute_reciprocal(T val, T n) {
    static_assert(sizeof(T) <= 8);
    return ((__uint128_t(1) << (n)) * val)/modulus;
}

template<typename T, T poly_dim, T modulus>
class NTT_engine {

    /* See for details
     * https://en.wikipedia.org/wiki/Root_of_unity_modulo_n#Finding_a_primitive_k-th_root_of_unity_modulo_n
https://crypto.stackexchange.com/questions/63614/finding-the-n-th-root-of-unity-in-a-finite-field 
*/
public:
    NTT_engine(): log_poly_dim(LOG_2<T, poly_dim>::n), r(LOG_2<T, modulus>::n), n(sizeof(T) * 8), m(r) {
        // is poly_dim a power of 2 ?
        static_assert((poly_dim & (poly_dim - 1)) == 0);
        // are there square roots for omega ?
        static_assert((modulus % (2 * poly_dim)) == 1);

        if (m <= n - 2) {
            s = r > 1 ? r - 2 : 0;
            t = n + 1;
        } else {
            s = r - 1;
            t = n;
        }

        // prime modulus carmichael = phi function
        // x^(carmichael(p)/k) = 1 ^ (-k) mod p
        //constexpr T cm = (modulus - 1) / poly_dim;
        constexpr T cm2 = (modulus - 1) / (2 * poly_dim);

        T tmp = 1;
        do {
            phi = mod_exp<T, modulus>(++tmp, cm2);

        } while(mod_exp<T, modulus>(phi, poly_dim) == 1);

        omega = mod_exp<T, modulus>(phi, 2);

        // use small fermat for inverse
        // a^(p-1) = 1 mod p => a^(p-2) = a^-1 mod p
        omega_inverse = mod_exp<T, modulus>(omega, modulus - 2);
        phi_inverse = mod_exp<T, modulus>(phi, modulus - 2);
        poly_dim_inverse = mod_exp<T, modulus>(poly_dim, modulus - 2);

        T o = 1, io = 1, p = 1, ip = 1;
        for(int i = 0; i < poly_dim; i++) {
            powers_of_omega[i] = o;
            powers_of_rec_omega[i] = compute_reciprocal<T, modulus>(o, n);
            powers_of_omega_inverse[i] = io;
            powers_of_rec_omega_inverse[i] = compute_reciprocal<T, modulus>(io, n);
            powers_of_phi[i] = p;
            powers_of_phi_inverse[i] = ip;


            o = mod_mul<T, modulus>(o, omega);
            io = mod_mul<T, modulus>(io, omega_inverse);
            p = mod_mul<T, modulus>(p, phi);
            ip = mod_mul<T, modulus>(ip, phi_inverse);
        }
    }

    void transform(std::array<T, poly_dim>& result, const std::array<T, poly_dim>& input) {
        for(int i = 0; i < poly_dim; i++) {
            forward_buffer[i] = mod_mul<T, modulus>(input[i], powers_of_phi[i]);
        }
        proto_transform(result, forward_buffer, powers_of_omega, powers_of_rec_omega);
    }

    void inverse_transform(std::array<T, poly_dim>& result, const std::array<T, poly_dim>& input) {
        proto_transform(forward_buffer, input, powers_of_omega_inverse, powers_of_rec_omega_inverse);

        for(int i = 0; i < poly_dim; i++) {
            result[i] = mod_mul<T, modulus>(forward_buffer[i], poly_dim_inverse);
            result[i] = mod_mul<T, modulus>(result[i], powers_of_phi_inverse[i]);
        }

    }

    void multiply(std::array<T, poly_dim>& result, const std::array<T, poly_dim>& lhs, const std::array<T, poly_dim>& rhs) {

        transform(product_buffer_1, rhs);
        transform(product_buffer_2, lhs);

        for(int i = 0; i < poly_dim; i++) {
            product_buffer_1[i] = mod_mul<T, modulus>(product_buffer_1[i], product_buffer_2[i]);
        }
        inverse_transform(result, product_buffer_1);
    }

private:
    // NOT PORTABLE
    void center_loop(T* result, T* buffer, uint32_t j, __m256i& twiddle, __m256i& rec_twiddle, const __m256i& qb) {

        //uint32_t print[8] = {0};

        const T* even_base = buffer + 2 * j;
        const __m256i permute_pattern = _mm256_setr_epi32(0,2,4,6,1,3,5,7);

        __m256i l1 = _mm256_loadu_si256((__m256i*)even_base);
        l1 = _mm256_permutevar8x32_epi32(l1, permute_pattern);
        __m256i l2 = _mm256_loadu_si256((__m256i*)(even_base + 8));
        l2 = _mm256_permutevar8x32_epi32(l2, permute_pattern);

        __m256i a = _mm256_permute2x128_si256(l1, l2, 0u | (2u << 4u));
        __m256i bc = _mm256_permute2x128_si256(l1, l2, 1u | (3u << 4u));

        /* compute b * twiddle mod Q, even slots*/
        __m256i b = _mm256_slli_epi64(bc, 32);
        b = _mm256_srli_epi64(b, 32);
        __m256i bv = _mm256_mul_epu32(b, rec_twiddle);
        bv = _mm256_srli_epi64(bv, 32);
        bv = _mm256_mul_epu32(bv, qb);
        __m256i bt = _mm256_mul_epi32(b, twiddle);
        bt = _mm256_sub_epi64(bt, bv);
        __m256i btmod = _mm256_min_epu32(_mm256_sub_epi32(bt, qb), bt);
        //

        /* compute b * twiddle mod Q, odd slots */
        b = _mm256_srli_epi64(bc, 32);
        bv = _mm256_mul_epu32(b, _mm256_srli_epi64(rec_twiddle, 32));
        bv = _mm256_srli_epi64(bv, 32);
        bv = _mm256_mul_epu32(bv, qb);
        bt = _mm256_mul_epi32(b, _mm256_srli_epi64(twiddle, 32));
        bt = _mm256_sub_epi64(bt, bv);
        __m256i btmod1 = _mm256_min_epu32(_mm256_sub_epi32(bt, qb), bt);
        //

        // merge btmod and btmod1
        btmod1 = _mm256_slli_si256(btmod1, 4);
        btmod = _mm256_xor_si256(btmod1, btmod);
        //////////////////////////////

        // compute sum
        __m256i apb = _mm256_add_epi32(a,btmod);
        apb = _mm256_min_epu32(apb, _mm256_sub_epi32(apb, qb));
        //////////////////////////////

        // compute diff
        __m256i cmp = _mm256_cmpgt_epi32(btmod, a);
        __m256i amb = _mm256_add_epi32(a, _mm256_and_si256(qb, cmp));
        amb = _mm256_sub_epi32(amb, btmod);
        //

        _mm256_storeu_si256((__m256i*)(result + j), apb);
        _mm256_storeu_si256((__m256i*)(result + poly_dim/2 + j), amb);
    }

    void proto_transform(std::array<T, poly_dim>& result, const std::array<T, poly_dim>& input, std::array<T, poly_dim>& twiddle_factors, std::array<T, poly_dim>& twiddle_factors_reciprocal) {


        for(uint32_t i = 0; i < poly_dim; i++) {
            result[i] = input[reverse_bits<uint32_t>(i) >> (sizeof(T) * 8 - log_poly_dim)];
        }

#ifdef OLD
        auto* half_result = transform_buffer.data() + poly_dim/2;
        for(uint32_t i = 0; i < log_poly_dim; i++) {
            uint32_t shift = log_poly_dim - i - 1;

            for(uint32_t j = 0; j < poly_dim/2; j++){

                uint32_t p = (j >> shift) << shift;
                T lhs = result[2 * j];
                T rhs = mod_mul<T, modulus>(result[2 * j + 1], twiddle_factors[p]);
                T inter = lhs + rhs;
                transform_buffer[j] = inter >= modulus ? inter - modulus : inter;
                half_result[j] = mod_sub<T, modulus>(lhs, rhs);
            }
            result = transform_buffer;
        }
#else
        __m256i twiddle1, twiddle2;// = _mm256_set1_epi32(powers_of_omega[0]);
        __m128i t_lo, t_hi;
        const __m256i qB = _mm256_set1_epi32(modulus);

        auto* t_buffer = result.data();
        auto* r_buffer = transform_buffer.data();

        for(uint32_t i = 0; i < log_poly_dim - 3 ; i++) {
            uint32_t shift = log_poly_dim - i - 1;
            for(uint32_t j = 0; j < poly_dim/2; j+= 8) {
                uint32_t p = (j >> shift) << shift;
                twiddle1 = _mm256_set1_epi32(twiddle_factors[p]);
                twiddle2 = _mm256_set1_epi32(twiddle_factors_reciprocal[p]);
                center_loop(r_buffer, t_buffer, j, twiddle1, twiddle2, qB);
            }
            std::swap(r_buffer, t_buffer);
        }

        for(uint32_t j = 0; j < poly_dim/2; j+=8) {
            t_lo = _mm_set1_epi32(twiddle_factors[j]);
            t_hi = _mm_set1_epi32(twiddle_factors[j + 4]);
            twiddle1 = _mm256_set_m128i(t_hi, t_lo);
            t_lo = _mm_set1_epi32(twiddle_factors_reciprocal[j]);
            t_hi = _mm_set1_epi32(twiddle_factors_reciprocal[j + 4]);
            twiddle2 = _mm256_set_m128i(t_hi, t_lo);
            center_loop(r_buffer, t_buffer, j, twiddle1, twiddle2, qB);
        }

        for(uint32_t j = 0; j < poly_dim/2; j+=8) {

            uint32_t lolo = twiddle_factors[j];
            uint32_t lohi = twiddle_factors[j + 2];
            uint32_t hilo = twiddle_factors[j + 4];
            uint32_t hihi = twiddle_factors[j + 6];

            twiddle1 = _mm256_setr_epi32(lolo, lolo, lohi, lohi, hilo, hilo, hihi, hihi);

            lolo = twiddle_factors_reciprocal[j];
            lohi = twiddle_factors_reciprocal[j + 2];
            hilo = twiddle_factors_reciprocal[j + 4];
            hihi = twiddle_factors_reciprocal[j + 6];

            twiddle2 = _mm256_setr_epi32(lolo, lolo, lohi, lohi, hilo, hilo, hihi, hihi);

            center_loop(t_buffer, r_buffer, j, twiddle1, twiddle2, qB);
        }

        for(uint32_t j = 0; j < poly_dim/2; j+=8) {
            twiddle1 = _mm256_loadu_si256((__m256i*)(twiddle_factors.data() + j));
            twiddle2 = _mm256_loadu_si256((__m256i*)(twiddle_factors_reciprocal.data() + j));
            center_loop(r_buffer, t_buffer, j, twiddle1, twiddle2, qB);
        }

#endif
    }

    T omega;
    T omega_inverse;
    T phi;
    T phi_inverse;
    T poly_dim_inverse;
    T log_poly_dim;

    /* barret params */
    T r, s, t, n, m, rec_p;

    alignas(32) std::array<T, poly_dim> transform_buffer;
    alignas(32) std::array<T, poly_dim> forward_buffer;
    alignas(32) std::array<T, poly_dim> product_buffer_1;
    alignas(32) std::array<T, poly_dim> product_buffer_2;

    std::array<T, poly_dim> powers_of_omega;
    alignas(16) std::array<T, poly_dim> powers_of_rec_omega;
    std::array<T, poly_dim> powers_of_omega_inverse;
    alignas(16) std::array<T, poly_dim> powers_of_rec_omega_inverse;
    std::array<T, poly_dim> powers_of_phi;
    std::array<T, poly_dim> powers_of_phi_inverse;
};

#endif