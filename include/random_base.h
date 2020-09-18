//
// Created by leonard on 10.09.20.
//

#ifndef LWE_STRUCTS_RANDOM_BASE_H
#define LWE_STRUCTS_RANDOM_BASE_H

#include <random>

template<typename T>
struct RandomEngine {
    virtual ~RandomEngine() = default;
    virtual T generate_normal() = 0;
    virtual T generate_uniform_lwe() = 0;
    virtual T generate_uniform_rlwe() = 0;
};

template<typename T>
struct DefaultRandomEngine : RandomEngine<T> {

    DefaultRandomEngine(T lwe_modulus, T ring_modulus, double normal_mean, double normal_stddev) : l_modulus(lwe_modulus), r_mod(ring_modulus), stddev(normal_stddev),
                                                                               mean(normal_mean), dev(), mt(dev()),
                                                                               normal_dist(normal_mean, normal_stddev),
                                                                               ring_uniform_dist(0, ring_modulus),
                                                                               lwe_uniform_dist(0, lwe_modulus) {};

    DefaultRandomEngine(const DefaultRandomEngine<T>& other) : DefaultRandomEngine(other.l_modulus, other.r_mod, other.mean, other.stddev) {

    }

    virtual T generate_normal()  {
        return normal_dist(mt);
    }

    virtual T generate_uniform_lwe() {
        return lwe_uniform_dist(mt);
    }

    virtual T generate_uniform_rlwe() {
        return ring_uniform_dist(mt);
    }


private:
    T l_modulus = 0;
    T r_mod = 0;
    double stddev = 0.0;
    double mean = 0.0;

    std::random_device dev;
    std::mt19937 mt;
    std::normal_distribution<double> ring_uniform_dist;
    std::normal_distribution<double> normal_dist;
    std::uniform_int_distribution<T> lwe_uniform_dist;
};

#endif //LWE_STRUCTS_RANDOM_BASE_H
