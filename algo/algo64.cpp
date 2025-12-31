/*
Copyright 2025, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

#include <gmp.h>

#define CHECK	true

inline void mpz_set_u64(mpz_t & rop, const uint64_t op)
{
	// _LONG_LONG_LIMB is required
	if (op == 0) mpz_set_ui(rop, 0);
	else { mpz_set_ui(rop, 1); rop->_mp_d[0] = op; }
}

#ifdef CHECK
// 2^p ?= 1 mod p
inline bool prp_slow(const uint64_t p)
{
	uint64_t r = 1, b = 2;
	for (uint64_t e = p - 1; e != 0; e /= 2)
	{
		if (e % 2 == 1) r = uint64_t((__uint128_t(r) * b) % p);
		b = uint64_t((__uint128_t(b) * b) % p);
	}
	return (r == 1);
}
#endif

typedef uint8_t		uint_8;
typedef int8_t		int_8;
typedef uint32_t	uint32;
typedef uint64_t	uint64;

// OpenCL functions
inline uint32 atomic_inc(uint32 * const p) { const uint32 t = *p; (*p)++; return t; }
inline uint64 mul_hi(const uint64 x, const uint64 y) { return uint64((x * __uint128_t(y)) >> 64); }

inline int ilog2(const uint64 x) { return 63 - __builtin_clzll(x); }
inline bool bittest(const uint64 x, const int b) { return ((x & ((uint64)(1) << b)) != 0); }

inline uint64 add_mod(const uint64 a, const uint64 b, const uint64 p) { return a + b - ((a >= p - b) ? p : 0); }
inline uint64 dup_mod(const uint64 a, const uint64 p) { return add_mod(a, a, p); }
inline uint64 sub_mod(const uint64 a, const uint64 b, const uint64 p) { return a - b + ((a < b) ? p : 0); }

// Montgomery modular multiplication
inline uint64 mul_mod(const uint64 a, const uint64 b, const uint64 p, const uint64 q)
{
	const uint64 ab_l = a * b, ab_h = mul_hi(a, b);
	return sub_mod(ab_h, mul_hi(ab_l * q, p), p);
}

inline uint64 sqr_mod(const uint64 a, const uint64 p, const uint64 q) { return mul_mod(a, a, p, q); }

inline uint64 mul_mod_const(const uint64 a, const uint64 b, const uint64 p, const uint64 bq)
{
	const uint64 ab_h = mul_hi(a, b);
	return sub_mod(ab_h, mul_hi(a * bq, p), p);
}

// Conversion out of Montgomery form
inline uint64 Montgomery2int(const uint64 r, const uint64 p, const uint64 q)
{
	const uint64 mp = mul_hi(r * q, p);
	return (mp != 0) ? p - mp : 0;
}

// p * p_inv = 1 (mod 2^64) (Newton's method)
inline uint64 invert(const uint64 p)
{
	uint64 p_inv = 2 - p;
	uint64 prev; do { prev = p_inv; p_inv *= 2 - p * p_inv; } while (p_inv != prev);
	return p_inv;
}

// a^e mod p, left-to-right algorithm
inline uint64 pow_mod(const uint64 a, const uint64 e, const uint64 p, const uint64 q)
{
	const uint64 aq = a * q;
	uint64 r = a;
	for (int b = ilog2(e) - 1; b >= 0; --b)
	{
		r = sqr_mod(r, p, q);
		if (bittest(e, b)) r = mul_mod_const(r, a, p, aq);
	}
	return r;
}

// 2^{(p - 1)/2} ?= +/-1 mod p
inline bool prp(const uint64 p, const uint64 q, const uint64 one)
{
	const uint64 e = (p - 1) / 2;
	int b = ilog2(e) - 1;
	uint64 r = dup_mod(one, p);	// 2 = 1 + 1
	r = dup_mod(r, p);			// 2 * 2 = 2 + 2
	if (bittest(e, b)) r = dup_mod(r, p);
	for (--b; b >= 0; --b)
	{
		r = sqr_mod(r, p, q);
		if (bittest(e, b)) r = dup_mod(r, p);
	}
	return ((r == one) || (r == p - one));
}

static void check(const uint64 k , const uint32 b, const int n)
{
	mpz_t zp, zr, zt; mpz_inits(zp, zr, zt, nullptr);

	mpz_set_u64(zp, k); mpz_mul_2exp(zp, zp, mp_bitcnt_t(n + 1)); mpz_add_ui(zp, zp, 1);
	char str[32]; mpz_get_str(str, 10, zp);

	mpz_set_ui(zr, b); mpz_powm_ui(zr, zr, uint32_t(1) << n, zp); mpz_add_ui(zr, zr, 1);
	if (mpz_cmp(zr, zp) == 0)
	{
		std::cout << str << " | " << b << "^{2^" << n << "}+1" << std::endl;
	}
	else
	{
		mpz_set_ui(zt, 2); mpz_powm(zt, zt, zp, zp);
		if (mpz_cmp_ui(zt, 2) != 0)
		{
			std::ostringstream ss; ss << str << " is not 2-prp." << std::endl;
			throw std::runtime_error(ss.str());
		}
		bool isPrime = true;
		for (uint32_t a = 3; a < 1000; a += 2)
		{
			mpz_set_ui(zt, a); mpz_powm(zt, zt, zp, zp);
			if (mpz_cmp_ui(zt, a) != 0)
			{
				isPrime = false;
				break;
			}
		}
		if (isPrime)
		{
			std::ostringstream ss; ss << str << " doesn't divide " << b << "^{2^" << n << "}+1." << std::endl;
			throw std::runtime_error(ss.str());
		}
	}

	mpz_clears(zp, zr, zt, nullptr);
}

static void test(const uint64 i_min, const uint64 i_max, const int n, const int log2_block_size, const uint_8 * const wheel, const int_8 * const kro_vector)
{
	const size_t block_size = size_t(1) << log2_block_size;
	const size_t factors_block = size_t(1) << 10;
	const size_t N_2_factors_block = (size_t(1) << (n - 1)) / factors_block;
	const int g_n = n + 1;

	uint64 * const k_vector = new uint64[block_size];
	uint64 * const q_vector = new uint64[block_size];
	uint64 * const ext_vector = new uint64[block_size];
	uint64 * const c_vector = new uint64[block_size];
	uint64 * const cn_vector = new uint64[block_size];

	uint32 prime_count = 0;

	for (uint64 i = i_min; i < i_max; ++i)
	{
		// generate_primes
		for (uint64 id = 0; id < block_size; ++id)
		{
			const uint64 j = (i << log2_block_size) | id;
			const uint64 k = 15 * (j / 8) + wheel[j % 8];

			// one is the Montgomery form of 1: 2^64 mod p = (2^64 - p) mod p
			const uint64 p = (k << g_n) | 1, q = invert(p), one = (-p) % p;
#ifdef CHECK
			if ((p % 3 == 0) || (p % 5 == 0)) throw std::runtime_error("Error: wheel.");
#endif
			const bool isprp = prp(p, q, one);
#ifdef CHECK
			if (isprp != prp_slow(p)) throw std::runtime_error("Error: prp.");
#endif
			if (isprp)
			{
				const uint32 prime_index = atomic_inc(&prime_count);
				k_vector[prime_index] = k;
				q_vector[prime_index] = q;
				ext_vector[prime_index] = one;
			}
		}

		// p = 1P: expected = 15/8 * 2 / log(1e15) * block_size = 0.1086 * block_size
		// p = 18446P: expected = 15/8 * 2 / log(18446e15) * block_size = 0.0845 * block_size
		const double pf = 15 / 8.0 * i * std::pow(2.0, double(log2_block_size + n + 1));
		const double expected = 15 / 8.0 * std::pow(2.0, log2_block_size + 1) / log(pf);
		std::cout << i << ": " << prime_count << " primes (" << expected << ")" << std::endl;

		// init_factors
		for (uint64 id = 0; id < block_size; ++id)
		{
			if (id >= prime_count) break;

			const uint64 k = k_vector[id], q = q_vector[id], one = ext_vector[id];

			const uint64 p = (k << g_n) | 1;
			const uint64 two = dup_mod(one, p);

			// p = 1 (mod 4). If a is odd then (a/p) = (p/a) = ({p mod a}/a)

			uint32 a = 3; uint64 am = add_mod(two, one, p);
			if (p % 3 != 2)
			{
				a += 2; am = add_mod(am, two, p);
				if (kro_vector[256 * ((5 - 3) / 2) + (uint32)(p % 5)] >= 0)
				{
					a += 2; am = add_mod(am, two, p);
					if (kro_vector[256 * ((7 - 3) / 2) + (uint32)(p % 7)] >= 0)
					{
						a += 4; am = add_mod(am, dup_mod(two, p), p);
						while (a < 256)
						{
							if (kro_vector[256 * ((a - 3) / 2) + (uint32)(p % a)] < 0) break;
							a += 2; am = add_mod(am, two, p);
						}
						if (a >= 256) { a = 0; am = 0; }
					}
				}
			}
#ifdef CHECK
			if (pow_mod(am, (p - 1) / 2, p, q) != p - one) throw std::runtime_error("Error: a.");
#endif
			// a^{(p - 1)/2} = -1 <=> a^{k*2^n} = -1. (a^k)^{2*i + 1} are the roots of b^{2^n} + 1 = 0 (mod p)
			const uint64 cm = pow_mod(am, k, p, q);
			const uint64 c = Montgomery2int(cm, p, q);
			c_vector[id] = c;
			ext_vector[id] = sqr_mod(cm, p, q);
			cn_vector[id] = p - c;
		}

		// check_factors
		for (size_t j = 0; j < N_2_factors_block; ++j)
		{
			// std::cout << j << " / " << N_2_factors_block << std::endl;

			for (uint64 id = 0; id < block_size; ++id)
			{
				if (id >= prime_count) break;

				uint64 c = c_vector[id];
				if (c == 0) continue;

				const uint64 k = k_vector[id], p = (k << g_n) | 1, q = q_vector[id];
				const uint64 c0sq = ext_vector[id], c0sqxq = c0sq * q;

#ifdef CHECK
				uint64 c0sq_check = Montgomery2int(c0sq, p, q);
#endif
				for (size_t l = 0; l < factors_block; ++l)
				{
					const uint64 b = (c % 2 == 0) ? c : p - c;
					const bool found = (b <= 2000000000);

					if (found) check(k, uint32(b), n);
#ifdef CHECK
					const uint64 cp = c;
#endif
					c = mul_mod_const(c, c0sq, p, c0sqxq);	// c = a^{(2*i + 1).k}
#ifdef CHECK
					if (c != uint64_t((__uint128_t(cp) * c0sq_check) % p)) throw std::runtime_error("Error: factors loop.");
#endif
				}

				if (j != N_2_factors_block - 1) c_vector[id] = c;
				else if (c != cn_vector[id]) throw std::runtime_error("Error: check failed.");
			}
		}

		prime_count = 0;
	}

	delete[] k_vector;
	delete[] q_vector;
	delete[] ext_vector;
	delete[] c_vector;
	delete[] cn_vector;
}

int main()
{
	try
	{
		// Kronecker symbols (i/j) for odd j <= 255.
		int_8 kro_vector[128 * 256];
		mpz_t zj; mpz_init(zj);
		for (uint32_t j = 3; j < 256; j += 2)
		{
			mpz_set_ui(zj, j);
			for (uint32_t i = 0; i < j; ++i)
			{
				kro_vector[256 * ((j - 3) / 2) + i] = int_8(mpz_ui_kronecker(i, zj));
			}
		}
		mpz_clear(zj);

		const int n = 23;

		// gcd(p mod 15, 15) = 1: 8 solutions
		uint_8 wheel[8];
		size_t i = 0;
		for (uint64 k = 0; k < 15; ++k)
		{
			const uint64 p = (k << (n + 1)) | 1;
			if ((p % 3 == 0) || (p % 5 == 0)) continue;
			wheel[i] = uint_8(k);
			++i;
		}
		if (i != 8) throw std::runtime_error("Error: wheel count.");

		const int log2_block_size = 12;
		const double unit = 1e15;

		const uint32 p_min = 100, p_max = p_min + 1;	// 1 <= p < 18446

		const double f = unit * 8.0 / 15 / std::pow(2.0, double(log2_block_size + n + 1));
		const uint64 i_min = uint64(std::floor(p_min * f)), i_max = uint64(std::ceil(p_max * f));

		mpz_t zp_min, zp_max; mpz_inits(zp_min, zp_max, nullptr);

		const uint64 j_min = i_min << log2_block_size, j_max = i_max << log2_block_size;
		const uint64 k_min = 15 * (j_min / 8) + wheel[j_min % 8], k_max = 15 * (j_max / 8) + wheel[j_max % 8];

		mpz_set_u64(zp_min, k_min); mpz_mul_2exp(zp_min, zp_min, n + 1); mpz_add_ui(zp_min, zp_min, 1);
		mpz_set_u64(zp_max, k_max); mpz_mul_2exp(zp_max, zp_max, n + 1); mpz_add_ui(zp_max, zp_max, 1);

		char p_min_str[32], p_max_str[32]; mpz_get_str(p_min_str, 10, zp_min); mpz_get_str(p_max_str, 10, zp_max);
		std::cout << "Testing n = " << n << ", p in [" << p_min_str << ", " << p_max_str << "[." << std::endl;

		mpz_clears(zp_min, zp_max, nullptr);

		std::cout << "For i = " << i_min << " to " << i_max - 1 << std::endl;

		test(i_min, i_max, n, log2_block_size, wheel, kro_vector);
	}
	catch (const std::runtime_error & e)
	{
		std::cout << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}