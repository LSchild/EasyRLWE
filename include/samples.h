//
// Created by leonard on 08.09.20.
//

#ifndef LWE_STRUCTS_SAMPLES_H
#define LWE_STRUCTS_SAMPLES_H

#include <cstddef>
#include "common.h"

template<typename T, size_t len, T modulus>
struct TemplateLweSample {

    TemplateLweSample() = default;

    TemplateLweSample(const TemplateVector<T, len, modulus>& oa, T ob) : a(oa), b(ob){ }

    TemplateLweSample(const TemplateLweSample& other) {
        a = other.a;
        b = other.b;
    }

    friend void swap(TemplateLweSample& lhs, TemplateLweSample& rhs) {
        using std::swap;

        swap(lhs.a, rhs.a);
        swap(lhs.b, rhs.b);
    }

    TemplateLweSample& operator=(TemplateLweSample other) {
        swap(*this, other);
        return *this;
    }

    TemplateVector<T, len, modulus> getA() const {
        return a;
    }

    T getB() const {
        return b;
    }

    T& operator[](int idx) {
        return a[idx];
    }

    T operator[](int idx) const {
        return a[idx];
    }

    void operator+=(const TemplateLweSample& rhs) {
        for(int i = 0; i < len; i++) {
            a[i] = (a[i] + rhs.a[i]) % modulus;
        }

        b = (b + rhs.b) % modulus;
    }

    void operator-=(const TemplateLweSample& rhs) {
        for(int i = 0; i < len; i++) {
            a[i] = mod_sub<T, modulus>(a[i], rhs.a[i]);
        }

        b = mod_sub<T, modulus>(b, rhs.b);
    }

    void operator*=(const T rhs) {
        for(int i = 0; i < len; i++) {
            a[i] = mod_mul<T, modulus>(a[i], rhs);
        }
        b = mod_mul<T, modulus>(b, rhs);
    }

    ~TemplateLweSample() = default;

private:
    TemplateVector<T, len, modulus> a{};
    T b = 0;
};

template<typename T, size_t dim, T modulus>
struct TemplateRLweSample {

    TemplateRLweSample() = default;

    TemplateRLweSample(const TemplateRLweSample& other) {
        a = other.a;
        b = other.b;
    }

    TemplateRLweSample(TemplatePolynomial<T, dim, modulus>& a0, TemplatePolynomial<T, dim, modulus>& a1) {
        a = a0;
        b = a1;
    }

    friend void swap(TemplateRLweSample& lhs, TemplateRLweSample& rhs) {
        using std::swap;

        swap(lhs.a, rhs.a);
        swap(lhs.b, rhs.b);
    }

    std::vector<TemplateRLweSample> decompose(T base) {
        std::vector<TemplateRLweSample> ret;

        auto lhs = a.decompose(base);
        auto rhs = b.decompose(base);

        for(int i = 0; i < lhs.size(); i++) {
            ret.emplace_back(lhs.at(i), rhs.at(i));
        }

        return ret;
    }

    TemplateRLweSample& operator=(TemplateRLweSample rhs) {
        swap(*this, rhs);

        return *this;
    }

    void operator+=(const TemplateRLweSample& rhs) {
        a += rhs.a;
        b += rhs.b;
    }

    void operator-=(const TemplateRLweSample& rhs) {
        a -= rhs.a;
        b -= rhs.b;
    }

    void operator*=(const TemplatePolynomial<T, dim, modulus>& rhs) {
        a *= rhs;
        b *= rhs;
    }

    TemplatePolynomial<T,dim,modulus> getA() const {
        return a;
    }

    TemplatePolynomial<T,dim,modulus> getB() const {
        return b;
    }

private:

    TemplatePolynomial<T, dim, modulus> a;
    TemplatePolynomial<T, dim, modulus> b;
};

template<typename T, size_t dim, T modulus>
struct TemplateRLwePrimeSample {

    TemplateRLwePrimeSample() = default;

    TemplateRLwePrimeSample(T dec_base, T dec_dim) : decomposition_base(dec_base), decomposition_dimension(dec_dim) {
        entries.reserve(dec_dim);
    }

    TemplateRLwePrimeSample(const TemplateRLwePrimeSample& other): entries(other.entries), decomposition_base(other.decomposition_base),
    decomposition_dimension(other.decomposition_dimension) {}

    friend void swap(TemplateRLwePrimeSample& lhs, TemplateRLwePrimeSample& rhs) {
        using std::swap;

        swap(lhs.entries, rhs.entries);
        swap(lhs.decomposition_dimension, rhs.decomposition_dimension);
        swap(lhs.decomposition_base, rhs.decomposition_base);

    }

    T getBase() const {
        return decomposition_base;
    }

    T getDimension() const {
        return decomposition_dimension;
    }

    TemplateRLwePrimeSample& operator=(TemplateRLwePrimeSample other) {
        swap(*this, other);
        return *this;
    }

    TemplateRLweSample<T, dim, modulus>& operator[](int idx) {
#ifdef DEBUG
        return entries.at(idx);
#else
        return entries[idx];
#endif
    }

    TemplateRLweSample<T, dim, modulus> operator[](int idx) const {
#ifdef DEBUG
        return entries.at(idx);
#else
        return entries[idx];
#endif
    }

    void operator+=(const TemplateRLwePrimeSample& rhs) {
        for(int i = 0; i < decomposition_dimension; i++) {
            entries.at(i) += rhs.entries.at(i);
        }
    }

    void operator-=(const TemplateRLwePrimeSample& rhs) {
        for(int i = 0; i < decomposition_dimension; i++) {
            entries[i] -= rhs.entries[i];
        }
    }

    std::vector<TemplateRLweSample<T, dim, modulus>>& getEntries() {
        return entries;
    }

private:
    std::vector<TemplateRLweSample<T, dim, modulus>> entries;
    T decomposition_base;
    T decomposition_dimension;
};

template<typename T, size_t dim, T modulus>
struct TemplateRGswSample {

    TemplateRGswSample() = default;

    TemplateRGswSample(const TemplateRGswSample& other): lhs(other.lhs), rhs(other.rhs) {
    }

    TemplateRGswSample(T dec_base, T dec_dim) : lhs(dec_base, dec_dim), rhs(dec_base, dec_dim) {
    }

    friend void swap(TemplateRGswSample& o_lhs, TemplateRGswSample& o_rhs) {
        using std::swap;

        swap(o_lhs.rhs, o_lhs.rhs);
        swap(o_lhs.lhs, o_rhs.lhs);

    }

    TemplateRGswSample& operator=(TemplateRGswSample other) {
        swap(*this, other);
        return *this;
    }

    void operator+=(const TemplateRGswSample& other) {
        lhs += other.lhs;
        rhs += other.rhs;
    }

    void operator-=(const TemplateRGswSample& other) {
        lhs -= other.lhs;
        rhs -= other.rhs;
    }

    TemplateRLwePrimeSample<T, dim, modulus>& operator[](int idx) {
        switch (idx) {
            case 0: return lhs;
            case 1: return rhs;
            default: throw std::invalid_argument("Indices greater than 1 not supported");
        }
    }

    TemplateRLwePrimeSample<T, dim, modulus> operator[](int idx) const {
        switch (idx) {
            case 0: return lhs;
            case 1: return rhs;
            default: throw std::invalid_argument("Indices greater than 1 not supported");
        }
    }

    TemplateRLwePrimeSample<T, dim, modulus>& getLHS() {
        return lhs;
    }

    TemplateRLwePrimeSample<T, dim, modulus>& getRHS() {
        return rhs;
    }

    TemplateRLwePrimeSample<T, dim, modulus> getLHS() const {
        return lhs;
    }

    TemplateRLwePrimeSample<T, dim, modulus> getRHS() const {
        return rhs;
    }


private:
    TemplateRLwePrimeSample<T, dim, modulus> lhs;
    TemplateRLwePrimeSample<T, dim, modulus> rhs;
};

template<typename T, size_t dim, T modulus>
TemplateLweSample<T, dim, modulus> operator+(TemplateLweSample<T, dim, modulus> lhs, const TemplateLweSample<T, dim, modulus>& rhs) {
    lhs += rhs;
    return lhs;
}

template<typename T, size_t dim, T modulus>
TemplateLweSample<T, dim, modulus> operator-(TemplateLweSample<T, dim, modulus> lhs, const TemplateLweSample<T, dim, modulus>& rhs) {
    lhs -= rhs;
    return lhs;
}

template<typename T, size_t dim, T modulus>
TemplateLweSample<T, dim, modulus> operator*(TemplateLweSample<T, dim, modulus> lhs, const T rhs) {
    lhs *= rhs;
    return lhs;
}

template<typename T, size_t dim, T modulus>
TemplateLweSample<T, dim, modulus> operator*(const T lsh, TemplateLweSample<T, dim, modulus> lhs) {
    lhs *= lsh;
    return lhs;
}

template<typename T, size_t dim, T modulus>
TemplateRLweSample<T, dim, modulus> operator+(TemplateRLweSample<T, dim, modulus> lhs, const TemplateRLweSample<T, dim, modulus>& rhs) {
    lhs += rhs;
    return lhs;
}

template<typename T, size_t dim, T modulus>
TemplateRLweSample<T, dim, modulus> operator-(TemplateRLweSample<T, dim, modulus> lhs, const TemplateRLweSample<T, dim, modulus>& rhs) {
    lhs -= rhs;
    return lhs;
}

template<typename T, size_t dim, T modulus>
TemplateRLweSample<T, dim, modulus> operator*(TemplateRLweSample<T, dim, modulus> lhs, const TemplatePolynomial<T, dim, modulus>& rhs) {
    lhs *= rhs;
    return lhs;
}

template<typename T, size_t dim, T modulus>
TemplateRLweSample<T, dim, modulus> operator*(const TemplatePolynomial<T, dim, modulus>& rhs, TemplateRLweSample<T, dim, modulus> lhs) {
    lhs *= rhs;
    return lhs;
}

template<typename T, size_t dim, T modulus>
TemplateRLwePrimeSample<T, dim, modulus> operator+(TemplateRLwePrimeSample<T, dim, modulus> lhs, const TemplateRLwePrimeSample<T, dim, modulus>& rhs) {
    lhs += rhs;
    return lhs;
}

template<typename T, size_t dim, T modulus>
TemplateRLwePrimeSample<T, dim, modulus> operator-(TemplateRLwePrimeSample<T, dim, modulus> lhs, const TemplateRLwePrimeSample<T, dim, modulus>& rhs) {
    lhs -= rhs;
    return lhs;
}

template<typename T, size_t dim, T modulus>
TemplateRLweSample<T, dim, modulus> operator*(TemplateRLwePrimeSample<T, dim, modulus> lhs, const TemplatePolynomial<T, dim, modulus>& rhs) {

    auto decomposed = rhs.decompose(lhs.getBase());
    TemplateRLweSample<T, dim, modulus> accu;

    for(int i = 0; i < lhs.getDimension(); i++) {
        accu += lhs[i] * decomposed[i];
    }

    return accu;
}

template<typename T, size_t dim, T modulus>
TemplateRLwePrimeSample<T, dim, modulus> operator*(const TemplatePolynomial<T, dim, modulus>& rhs, TemplateRLwePrimeSample<T, dim, modulus> lhs) {
    return lhs * rhs;
}

template<typename T, size_t dim, T modulus>
TemplateRGswSample<T, dim, modulus> operator+(TemplateRGswSample<T, dim, modulus> lhs, const TemplateRGswSample<T, dim, modulus>& rhs) {
    lhs += rhs;
    return lhs;
}

template<typename T, size_t dim, T modulus>
TemplateRGswSample<T, dim, modulus> operator-(TemplateRGswSample<T, dim, modulus> lhs, const TemplateRGswSample<T, dim, modulus>& rhs) {
    lhs -= rhs;
    return lhs;
}

template<typename T, size_t dim, T modulus>
TemplateRLweSample<T, dim, modulus> operator*(TemplateRGswSample<T, dim, modulus> lhs, const TemplateRLweSample<T, dim, modulus>& rhs) {
    return (lhs.getLHS() * rhs.getA()) + (lhs.getRHS() * rhs.getB());
}

template<typename T, size_t dim, T modulus>
TemplateRLweSample<T, dim, modulus> operator*( const TemplateRLweSample<T, dim, modulus>& rhs, TemplateRGswSample<T, dim, modulus> lhs ) {
    return lhs * rhs;
}

#endif //LWE_STRUCTS_SAMPLES_H
