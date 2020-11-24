//
// Created by leonard on 13.11.20.
//

#ifndef LWE_STRUCTS_VECTOR_OPS_H
#define LWE_STRUCTS_VECTOR_OPS_H

#include <cstddef>
#include <immintrin.h>

// assumes AVX 2 support, modulus < 2^type_bitwidth - 1

template<typename T, size_t dim, T modulus>
void modadd_vec(T* result, T* lhs, const T* rhs) {

    static_assert(sizeof(T) <= 8);

    if constexpr(sizeof(T) == 4) {

        __m256i q8 = _mm256_set1_epi32(modulus);
        __m256i t8 = _mm256_sub_epi32(q8, _mm256_set1_epi32(1));

        for(int i = 0; i < dim; i+=8) {
            __m256i a8 = _mm256_loadu_si256((__m256i*)&lhs[i]);
            __m256i b8 = _mm256_loadu_si256((__m256i*)&rhs[i]);
            __m256i r8 = _mm256_add_epi32(a8,b8);
            __m256i cmp = _mm256_cmpgt_epi32(r8, t8);

            // cmp = _mm256_and_si256(r8, cmp);
            r8 = _mm256_sub_epi32(r8, _mm256_and_si256(q8, cmp));
            _mm256_storeu_si256((__m256i*)&result[i], r8);
        }

    }  else {

        __m256i q8 = _mm256_set1_epi64x(modulus);
        __m256i t8 = _mm256_sub_epi64(q8, _mm256_set1_epi32(1));

        for(int i = 0; i < dim; i+=4) {
            __m256i a8 = _mm256_loadu_si256((__m256i*)&lhs[i]);
            __m256i b8 = _mm256_loadu_si256((__m256i*)&rhs[i]);


            __m256i r8 = _mm256_add_epi64(a8,b8);
            __m256i cmp = _mm256_cmpgt_epi64(r8, t8);

            // cmp = _mm256_and_si256(r8, cmp);
            r8 = _mm256_sub_epi64(r8, _mm256_and_si256(q8, cmp));


            _mm256_storeu_si256((__m256i*)&result[i], r8);
        }
    }
}

template<typename T, size_t dim, T modulus>
void modsub_vec(T* result, T* lhs, const T* rhs) {

    static_assert(sizeof(T) <= 8);
    static_assert(dim % 8 == 0);

    if constexpr (sizeof(T) == 4) {

        __m256i q8 = _mm256_set1_epi32(modulus);

        for(int i = 0; i < dim; i+=8) {
            __m256i a8 = _mm256_loadu_si256((__m256i*)&lhs[i]);
            __m256i b8 = _mm256_loadu_si256((__m256i*)&rhs[i]);
            __m256i cmp = _mm256_cmpgt_epi32(b8, a8);
            __m256i s = _mm256_add_epi32(a8, _mm256_and_si256(q8, cmp));
            __m256 r = _mm256_sub_epi32(s, b8);
            _mm256_storeu_si256((__m256i*)&result[i], r);
        }

    } else {

        __m256i q8 = _mm256_set1_epi64x(modulus);

        for(int i = 0; i < dim; i+=8) {
            __m256i a8 = _mm256_loadu_si256((__m256i*)&lhs[i]);
            __m256i b8 = _mm256_loadu_si256((__m256i*)&rhs[i]);
            __m256i cmp = _mm256_cmpgt_epi64(b8, a8);
            __m256i s = _mm256_add_epi64(a8, _mm256_and_si256(q8, cmp));
            __m256 r = _mm256_sub_epi64(s, b8);
            _mm256_storeu_si256((__m256i*)&result[i], r);
        }
    }

}

#endif //LWE_STRUCTS_VECTOR_OPS_H
