#ifndef COMMON_HPP
#define COMMON_HPP

#ifndef __CUDA_ARCH__
#define __device__
#define __host__
#endif

#define __CLANG__

template<typename T, T modulus>
__device__ __host__ inline T mod_sub(T a, T b) {
	return (a >= b) ? a - b : (modulus - b) + a;
}

// maybe inline
template<typename T, T modulus>
__device__ __host__ T mod_mul(T a, T b) {

	if constexpr (sizeof(T) <= 4) {
		return (uint64_t(a) * uint64_t(b)) % modulus;
	} else if constexpr (sizeof(T) == 8) {
		return (__uint128_t(a) * __uint128_t(b)) % modulus;
	} else {
		T res = 0;
		while (b > 0) {
			res = (res + a * (b & 1)) % modulus;
			b >>= 1;
			a = (2 * a) % modulus;
		}
		return res;
	}
}

//maybe inline
template<typename T, T modulus>
__device__ __host__ T mod_exp(T a, T p) {
	T res = 1;
	T accu = a;
	while (p > 0) {
		if (p & 1) {
			res = mod_mul<T, modulus>(res, accu);
		}
		accu = mod_mul<T, modulus>(accu, accu);
		p >>= 1;
	}
	return res;
}

// https://stackoverflow.com/questions/63776/bit-reversal-of-an-integer-ignoring-integer-size-and-endianness
template<typename T>
__device__ __host__ T reverse_bits_no_specials(T in) {
	T nrev = 0, i, bit1, bit2;

	int count;
	constexpr int type_bits = sizeof(T) * 8;

	for(i = 0; i < type_bits; i += 2)
	{
		/*In each iteration, we  swap one bit on the 'right half'
		  of the number with another on the left half*/

		count =  type_bits - i - 1;  /*this is used to find how many positions
					       to the left (and right) we gotta move
					       the bits in this iteration*/

		bit1 = in & (1<<(i/2)); /*Extract 'right half' bit*/
		bit1 <<= count;         /*Shift it to where it belongs*/

		bit2 = in & 1<<((i/2) + count);  /*Find the 'left half' bit*/
		bit2 >>= count;         /*Place that bit in bit1's original position*/

		nrev |= bit1;   /*Now add the bits to the reversal result*/
		nrev |= bit2;
	}
	return nrev;
}

template<typename T>
inline __device__ __host__ T reverse_bits(T in) {

#if defined(__CUDA_ARCH__)

	if constexpr(sizeof(T) == 4)
		return __brev(in);
	else if constexpr(sizeof(T) == 8)
		return __brevll(in);
	else
		return reverse_bits_no_specials<T>(in);

#elif defined(__CLANG__)
	if constexpr(sizeof(T) == 1)
		return __builtin_bitreverse8(in);
	else if constexpr(sizeof(T) == 2)
		return __builtin_bitreverse16(in);
	else if constexpr(sizeof(T) == 4)
		return __builtin_bitreverse32(in);
	else if constexpr(sizeof(T) == 8)
		return __builtin_bitreverse64(in);
	// maybe for future ints ? __uint128 is a thing after all...
	else
		return reverse_bits_no_specials<T>(in);
#else
	return reverse_bits_no_specials<T>(in);
#endif

}

template<typename T, T val>
struct PLOG_2 {
	const static T n = PLOG_2<T, (val >> 1)>::n + 1;
};

template<>
struct PLOG_2<uint32_t, 0> {
	const static uint32_t n = 0;
};

template<>
struct PLOG_2<uint64_t, 0> {
	const static uint64_t n = 0;
};

template<>
struct PLOG_2<__uint128_t, 0> {
	const static __uint128_t n = 0;
};

template<typename T, T val>
struct LOG_2 {
	const static T n = PLOG_2<T, val>::n - 1;
};

#endif
