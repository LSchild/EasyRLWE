//
// Created by leonard on 08.09.20.
//


#ifndef LWE_STRUCTS_COMMON_H
#define LWE_STRUCTS_COMMON_H

#include <cstddef>
#include <array>
#include <vector>
#include <cmath>

#include "ntt.h"

NTT_engine<uint32_t, 1024, 134215681> default_uint32_t_engine{};

enum Format {
    NTT,
    DEFAULT
};

template<typename T, size_t len, T mod>
struct TemplateVector {

    TemplateVector() = default;

    explicit TemplateVector(T a) {
        elements.fill(a);
    }

    TemplateVector(TemplateVector& other) : elements(other.elements) {};

    friend void swap(TemplateVector& lhs, TemplateVector& rhs) {
        using std::swap;
        swap(lhs.elements, rhs.elements);
    }

    TemplateVector& operator=(TemplateVector other) {
        swap(*this, other);
        return *this;
    }

    T& operator[](int idx) {
        return elements[idx];
    }

    T operator[](int idx) const {
        return elements[idx];
    }

    void operator+=(TemplateVector& rhs) {
        for(int i = 0; i < len; i++) {
            elements[i] = (elements[i] + rhs[i]) % mod;
        }
    }

    void operator-=(TemplateVector& rhs) {
        for(int i = 0; i < len; i++) {
            elements[i] = mod_sub<T, mod>(elements[i], rhs[i]);
        }
    }

    void operator*=(T scalar) {
        for(auto& v: elements) {
            v = mod_mul<T, mod>(v, scalar);
        }
    }

    const std::array<T, len>& data() {
        return elements;
    }

private:
    std::array<T, len> elements{};
};

template<typename T, size_t len, T mod>
TemplateVector<T, len, mod> operator+(TemplateVector<T, len, mod> lhs, const TemplateVector<T, len, mod>& rhs) {
    lhs += rhs;
    return lhs;
}
template<typename T, size_t len, T mod>
TemplateVector<T, len, mod> operator-(TemplateVector<T, len, mod> lhs, const TemplateVector<T, len, mod>& rhs) {
    lhs -= rhs;
    return lhs;
}

template<typename T, size_t len, T mod>
TemplateVector<T, len, mod> operator*(TemplateVector<T, len, mod> lhs, const T& rhs) {
    lhs *= rhs;
    return lhs;
}

/*
 * Dot product
 */
template<typename T, size_t len, T mod>
T operator*(const TemplateVector<T, len, mod>& lhs, const TemplateVector<T, len, mod>& rhs) {
    T accu(0);
    T temp(0);

    for(int i = 0; i < len; i++) {
        temp = mod_mul<T, mod>(lhs[i], rhs[i]);
        accu = (accu + temp) % mod;
    }
}

template<typename T, size_t dim, T mod>
struct TemplatePolynomial {

    TemplatePolynomial(): ntt(default_uint32_t_engine) {

    }

    explicit TemplatePolynomial(NTT_engine<T, dim, mod>& nttEngine) : ntt(nttEngine) {}

    TemplatePolynomial(const TemplatePolynomial& other) : ntt(other.ntt)  {
        this->coefficients = other.coefficients;
    }

    TemplatePolynomial(NTT_engine<T, dim, mod>& nttEngine, T val, int idx) : ntt(nttEngine) {
        coefficients[idx] = val;
    }

    explicit TemplatePolynomial(NTT_engine<T, dim, mod>& nttEngine, TemplateVector<T, dim, mod>& vector): ntt(nttEngine) {
        this->coefficients = vector.data();
    }

    friend void swap(TemplatePolynomial& lhs, TemplatePolynomial& rhs) {
        using std::swap;

        swap(lhs.coefficients, rhs.coefficients);
        swap(lhs.ntt, rhs.ntt);
    }

    std::vector<TemplateVector<T, dim, mod>> decompose(T base) const {

        std::vector<TemplateVector<T, dim, mod>> decomposed;
        TemplateVector copy(this);

        auto slots = (size_t) (std::ceil(std::logb(mod) / std::logb(base)));
        decomposed.reserve(slots);

        for(int i = 0; i < slots; i++) {
            decomposed.emplace_back(this->ntt);
        }

        if (base && ((base & (base - 1)) == 0)) {


            T mask = base - 1;
            T shift = T(std::logb(base));

            for(int i = 0; i < slots; i++) {
                for(int j = 0; j < slots; j++) {
                    decomposed[i][j] = copy[j] & mask;
                    copy[j] >> shift;
                }
            }

        } else {

            for(int i = 0; i < slots; i++) {
                for(int j = 0; j < slots; j++) {
                    decomposed[i][j] = copy[j] % base;
                    copy[j] /= base;
                }
            }

        }

        return decomposed;
    }

    void setFormat(Format target) {
        if(target == fmt) return;

        if(target == DEFAULT) {
            ntt.inverse_transform(coefficients, coefficients);
        } else {
            ntt.transform(coefficients, coefficients);
        }
    }

    TemplatePolynomial& operator=(TemplatePolynomial other) {
        swap(*this, other);
        return *this;
    }

    void operator+=(const TemplatePolynomial& other) {
        for(int i = 0; i < dim; i++) {
            coefficients[i] = (coefficients[i] + other.coefficients[i]) % mod;
        }
    }

    void operator-=(const TemplatePolynomial& other) {
        for(int i = 0; i < dim; i++) {
            T a = coefficients[i];
            T b = other.coefficients[i];

            coefficients[i] = (a >= b) ? a - b : (mod - b) + a;
        }
    }

    void operator*=(T rhs) {
        for(int i = 0; i < dim; i++) {
            coefficients[i] = mod_mul<T, mod>(coefficients[i], rhs);
        }
    }

    void operator*=(const TemplatePolynomial& rhs) {
#ifdef DEBUG
        if (rhs.fmt != this->fmt)
            throw std::invalid_argument("Format of this and rhs do not match !");
#endif
        if(rhs.fmt == DEFAULT)
            ntt.multiply(coefficients, coefficients, rhs.coefficients);
        else {
            for(int i = 0; i < dim; i++) {
                coefficients[i] = mod_mul<T, mod>(coefficients[i], rhs.coefficients[i]);
            }
        }
    }

    T& operator[](int idx) {
        return coefficients[idx];
    }

    T operator[](int idx) const {
        return coefficients[idx];
    }

private:

    std::array<T, dim> coefficients{};
    NTT_engine<T, dim, mod>& ntt;
    Format fmt = DEFAULT;

};

template<typename T,  size_t dim, T modulus>
TemplatePolynomial<T, dim, modulus> operator+(TemplatePolynomial<T, dim, modulus> lhs, const TemplatePolynomial<T, dim, modulus>& rhs) {
    lhs += rhs;
    return lhs;
}

template<typename T,  size_t dim, T modulus>
TemplatePolynomial<T, dim, modulus> operator*(TemplatePolynomial<T, dim, modulus> lhs, TemplatePolynomial<T, dim, modulus>& rhs) {
    lhs *= rhs;
    return lhs;
}

template<typename T,  size_t dim, T modulus>
TemplatePolynomial<T, dim, modulus> operator-(TemplatePolynomial<T, dim, modulus> lhs, const TemplatePolynomial<T, dim, modulus>& rhs) {
    lhs -= rhs;
    return lhs;
}

#endif //LWE_STRUCTS_COMMON_H
