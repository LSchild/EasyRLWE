//
// Created by leonard on 16.09.20.
//

#ifndef LWE_STRUCTS_PARAMETERS_H
#define LWE_STRUCTS_PARAMETERS_H

#include "lwe.h"

#define RLWE_DIM 1024
#define RLWE_MOD 134215681

#define LWE_DIM 512
#define LWE_MOD 1024

#define PT_DIM 32
#define DELTA (uint32_t(RLWE_MOD) / PT_DIM)

#define NORMAL_STDDEV 3.19
#define ELEMENT_TYPE uint32_t

using Vector = TemplateVector<ELEMENT_TYPE, LWE_DIM, LWE_MOD>;
using RingVector = TemplateVector<ELEMENT_TYPE, RLWE_DIM, RLWE_MOD>;
using Polynomial = TemplatePolynomial<ELEMENT_TYPE, RLWE_DIM, RLWE_MOD>;
using CryptoEngine = TemplateCryptoEngine<ELEMENT_TYPE, LWE_DIM, LWE_MOD, RLWE_DIM, RLWE_MOD>;

#endif //LWE_STRUCTS_PARAMETERS_H
