/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#ifdef __NV_CL_C_VERSION
	#define PTX_ASM	1
#endif

#if !defined(G_N)
#define G_N			16
#define LN_BLKSZ	22
#define FBLK		1024
#define uint_8		uchar

__constant uint_8 wheel[8] = { 0, 1, 3, 6, 7, 10, 12, 13 };
#endif

typedef uint	sz_t;
typedef char	int_8;
typedef uint	uint_32;
typedef ulong	uint_64;
typedef ulong2	uint_64_2;

inline int ilog2(const uint_64 x) { return 63 - clz(x); }
inline bool bittest(const uint_64 x, const int b) { return ((x & ((uint_64)(1) << b)) != 0); }

inline uint_64 add_mod(const uint_64 a, const uint_64 b, const uint_64 p) { return a + b - ((a >= p - b) ? p : 0); }
inline uint_64 dup_mod(const uint_64 a, const uint_64 p) { return add_mod(a, a, p); }
inline uint_64 sub_mod(const uint_64 a, const uint_64 b, const uint_64 p) { return a - b + ((a < b) ? p : 0); }

// Montgomery modular multiplication
inline uint_64 mul_mod(const uint_64 a, const uint_64 b, const uint_64 p, const uint_64 q)
{
	const uint_64 ab_l = a * b, ab_h = mul_hi(a, b);
	return sub_mod(ab_h, mul_hi(ab_l * q, p), p);
}

inline uint_64 sqr_mod(const uint_64 a, const uint_64 p, const uint_64 q) { return mul_mod(a, a, p, q); }

inline uint_64 mul_mod_const(const uint_64 a, const uint_64 b, const uint_64 p, const uint_64 bq)
{
	const uint_64 ab_h = mul_hi(a, b);
	return sub_mod(ab_h, mul_hi(a * bq, p), p);
}

// Conversion out of Montgomery form
inline uint_64 Montgomery2int(const uint_64 r, const uint_64 p, const uint_64 q)
{
	const uint_64 mp = mul_hi(r * q, p);
	return (mp != 0) ? p - mp : 0;
}

// p * p_inv = 1 (mod 2^64) (Newton's method)
inline uint_64 invert(const uint_64 p)
{
	uint_64 p_inv = 2 - p;
	uint_64 prev; do { prev = p_inv; p_inv *= 2 - p * p_inv; } while (p_inv != prev);
	return p_inv;
}

// a^e mod p, left-to-right algorithm
inline uint_64 pow_mod(const uint_64 a, const uint_64 e, const uint_64 p, const uint_64 q)
{
	const uint_64 aq = a * q;
	uint_64 r = a;
	for (int b = ilog2(e) - 1; b >= 0; --b)
	{
		r = sqr_mod(r, p, q);
		if (bittest(e, b)) r = mul_mod_const(r, a, p, aq);
	}
	return r;
}

// 2^{(p - 1)/2} ?= +/-1 mod p
inline bool prp(const uint_64 p, const uint_64 q, const uint_64 one)
{
	const uint_64 e = (p - 1) / 2;
	int b = ilog2(e) - 1;
	uint_64 r = dup_mod(one, p);	// 2 = 1 + 1
	r = dup_mod(r, p);			// 2 * 2 = 2 + 2
	if (bittest(e, b)) r = dup_mod(r, p);
	for (--b; b >= 0; --b)
	{
		r = sqr_mod(r, p, q);
		if (bittest(e, b)) r = dup_mod(r, p);
	}
	return ((r == one) || (r == p - one));
}

__kernel
void generate_primes(__global sz_t * restrict const prime_count,
	__global uint_64 * restrict const k_vector, __global uint_64 * restrict const q_vector,
	__global uint_64 * restrict const ext_vector, const uint_64 i)
{
	const sz_t id = (sz_t)get_global_id(0);

	const uint_64 j = (i << LN_BLKSZ) | id;
	const uint_64 k = 15 * (j / 8) + wheel[j % 8];

	// one is the Montgomery form of 1: 2^64 mod p = (2^64 - p) mod p
	const uint_64 p = (k << G_N) | 1, q = invert(p), one = (-p) % p;
	if (prp(p, q, one))
	{
		const uint prime_index = atomic_inc(prime_count);
		k_vector[prime_index] = k;
		q_vector[prime_index] = q;
		ext_vector[prime_index] = one;
	}
}

__kernel
void init_factors(__global const sz_t * restrict const prime_count,
	__global const uint_64 * restrict const k_vector, __global const uint_64 * restrict const q_vector,
	__global uint_64 * restrict const ext_vector, __global const int_8 * restrict const kro_vector,
	__global uint_64 * restrict const c_vector, __global uint_64 * restrict const cn_vector)
{
	const sz_t id = (sz_t)get_global_id(0);
	if (id >= *prime_count) return;

	const uint_64 k = k_vector[id], q = q_vector[id], one = ext_vector[id];

	const uint_64 p = (k << G_N) | 1;
	const uint_64 two = dup_mod(one, p);

	// p = 1 (mod 4). If a is odd then (a/p) = (p/a) = ({p mod a}/a)

	uint_32 a = 3; uint_64 am = add_mod(two, one, p);
	if (p % 3 != 2)
	{
		a += 2; am = add_mod(am, two, p);
		if (kro_vector[256 * ((5 - 3) / 2) + (uint_32)(p % 5)] >= 0)
		{
			a += 2; am = add_mod(am, two, p);
			if (kro_vector[256 * ((7 - 3) / 2) + (uint_32)(p % 7)] >= 0)
			{
				a += 4; am = add_mod(am, dup_mod(two, p), p);
				while (a < 256)
				{
					if (kro_vector[256 * ((a - 3) / 2) + (uint_32)(p % a)] < 0) break;
					a += 2; am = add_mod(am, two, p);
				}
				if (a >= 256)
				{
					c_vector[id] = 0;
					return;
				}
			}
		}
	}

	// a^{(p - 1)/2} = -1 <=> a^{k*2^n} = -1. (a^k)^{2*i + 1} are the roots of b^{2^n} + 1 = 0 (mod p)
	const uint_64 cm = pow_mod(am, k, p, q);
	const uint_64 c = Montgomery2int(cm, p, q);
	c_vector[id] = c;
	ext_vector[id] = sqr_mod(cm, p, q);
	cn_vector[id] = p - c;
}

__kernel
void check_factors(__global const sz_t * restrict const prime_count,
	__global const uint_64 * restrict const k_vector, __global const uint_64 * restrict const q_vector,
	__global uint_64 * restrict const c_vector, __global const uint_64 * restrict const ext_vector,
	__global const uint_64 * restrict const cn_vector,
	__global sz_t * restrict const factor_count, __global uint_64_2 * restrict const factor_vector,
	__global sz_t * restrict const error_count, __global uint_64 * restrict const error_vector,
	const char last)
{
	const sz_t id = (sz_t)get_global_id(0);
	if (id >= *prime_count) return;

	uint_64 c = c_vector[id];
	if (c == 0) return;

	const uint_64 k = k_vector[id], p = (k << G_N) | 1, q = q_vector[id];
	const uint_64 c0sq = ext_vector[id], c0sqxq = c0sq * q;

	for (size_t l = 0; l < FBLK; ++l)
	{
		const uint_64 b = (c % 2 == 0) ? c : p - c;

		if (b <= 2000000000)
		{
			const sz_t factor_index = atomic_inc(factor_count);
			factor_vector[factor_index] = (uint_64_2)(k, b);
		}

		c = mul_mod_const(c, c0sq, p, c0sqxq);	// c = a^{(2*i + 1).k}
	}

	if (last == (char)(0)) c_vector[id] = c;
	else if (c != cn_vector[id])
	{
		const sz_t error_index = atomic_inc(error_count);
		error_vector[error_index] = k;
	}
}

__kernel
void clear(__global sz_t * restrict const count)
{
	*count = 0;
}
