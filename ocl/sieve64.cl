/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

typedef uint	uint32;
typedef ulong	uint64;

inline int ilog2(const uint64 x) { return 63 - clz(x); }

inline bool bittest(const uint64 x, const int b) { return ((x & ((uint64)(1) << b)) != 0); }

inline uint64 add_mod(const uint64 a, const uint64 b, const uint64 p)
{
	const uint64 c = (a >= p - b) ? p : 0;
	return a + b - c;
}

inline uint64 sub_mod(const uint64 a, const uint64 b, const uint64 p)
{
	const uint64 c = (a < b) ? p : 0;
	return a - b + c;
}

inline uint64 mul_mod(const uint64 a, const uint64 b, const uint64 p, const uint64 q)
{
	const uint64 ab_l = a * b, ab_h = mul_hi(a, b);
	const uint64 mp = mul_hi(ab_l * q, p), r = ab_h - mp;
	return (ab_h < mp) ? r + p : r;
}

inline uint64 mul_mod_const(const uint64 a, const uint64 b, const uint64 p, const uint64 q, const uint64 bq)
{
	const uint64 ab_h = mul_hi(a, b);
	const uint64 mp = mul_hi(a * bq, p), r = ab_h - mp;
	return (ab_h < mp) ? r + p : r;
}

inline uint64 toInt(const uint64 r, const uint64 p, const uint64 q)
{
	const uint64 mp = mul_hi(r * q, p);
	return (mp != 0) ? p - mp : 0;
}

// p * p_inv = 1 (mod 2^64) (Newton's method)
inline uint64 invert(const uint64 p)
{
	uint64 p_inv = 1, prev = 0;
	while (p_inv != prev) { prev = p_inv; p_inv *= 2 - p * p_inv; }
	return p_inv;
}

// a^e mod p, left-to-right algorithm
inline uint64 pow_mod(const uint64 a, const uint64 e, const uint64 p, const uint64 q)
{
	const uint64 aq = a * q;
	uint64 r = a;
	for (int b = ilog2(e) - 1; b >= 0; --b)
	{
		r = mul_mod(r, r, p, q);
		if (bittest(e, b)) r = mul_mod_const(r, a, p, q, aq);
	}
	return r;
}

// 2^{(p - 1)/2} ?= +/-1 mod p
inline bool prp(const uint64 p, const uint64 q, const uint64 one)
{
	const uint64 e = (p - 1) / 2;
	int b = ilog2(e) - 1;
	uint64 r = add_mod(one, one, p); r = add_mod(r, r, p);
	if (bittest(e, b)) r = add_mod(r, r, p);
	for (--b; b >= 0; --b)
	{
		r = mul_mod(r, r, p, q);
		if (bittest(e, b)) r = add_mod(r, r, p);
	}
	return ((r == one) || (r == p - one));
}

__kernel
void generate_primes(__global uint * restrict const prime_count, __global ulong2 * restrict const p_vector,
	__global ulong2 * restrict const q_vector, __global ulong2 * restrict const one_vector, const ulong i)
{
	const uint64 k = (i << log2GlobalWorkSize) | get_global_id(0);

	const uint64 p = (k << (gfn_n + 1)) | 1, q = invert(p), one = (-p) % p;
	if (prp(p, q, one))
	{
		const uint prime_index = atomic_inc(prime_count);
		p_vector[prime_index] = (ulong2)(p, (ulong)(0));
		q_vector[prime_index] = (ulong2)(q, (ulong)(0));
		one_vector[prime_index] = (ulong2)(one, (ulong)(0));
	}
}

__kernel
void init_factors(__global const uint * restrict const prime_count, __global const ulong2 * restrict const p_vector,
	__global const ulong2 * restrict const q_vector, __global const ulong2 * restrict const one_vector,
	__global const char * restrict const _kro_vector,
	__global ulong2 * restrict const c_vector, __global ulong2 * restrict const a2k_vector)
{
	const size_t i = get_global_id(0);
	if (i >= *prime_count) return;

	const uint64 p = p_vector[i].s0, q = q_vector[i].s0, one = one_vector[i].s0;
	const uint64 two = add_mod(one, one, p);

	// p = 1 (mod 4). If a is odd then (a/p) = (p/a) = ({p mod a}/a)

	uint32 a = 3; uint64 am = add_mod(two, one, p);
	if (p % 3 != 2)
	{
		a += 2; am = add_mod(am, two, p);
		if (_kro_vector[256 * ((5 - 3) / 2) + (uint32)(p % 5)] >= 0)
		{
			a += 2; am = add_mod(am, two, p);
			if (_kro_vector[256 * ((7 - 3) / 2) + (uint32)(p % 7)] >= 0)
			{
				a += 4; am = add_mod(am, two, p); am = add_mod(am, two, p);
				while (a < 256)
				{
					if (_kro_vector[256 * ((a - 3) / 2) + (uint32)(p % a)] < 0) break;
					a += 2; am = add_mod(am, two, p);
				}
				if (a >= 256)
				{
					c_vector[i] = (ulong2)(0, 0);
					return;
				}
			}
		}
	}

	const uint64 k = p >> (gfn_n + 1);
	const uint64 cm = pow_mod(am, k, p, q);

	const uint64 c = toInt(cm, p, q);
	c_vector[i] = (ulong2)(c, 0);

	const uint64 a2km = mul_mod(cm, cm, p, q);
	a2k_vector[i] = (ulong2)(a2km, 0);
}

__kernel
void check_factors(__global const uint * restrict const prime_count,
	__global const ulong2 * restrict const p_vector, __global const ulong2 * restrict const q_vector,
	__global ulong2 * restrict const c_vector, __global const ulong2 * restrict const a2k_vector,
	__global uint * restrict const factor_count, __global ulong2 * restrict const factor)
{
	const size_t i = get_global_id(0);
	if (i >= *prime_count) return;

	const uint64 p = p_vector[i].s0, q = q_vector[i].s0;

	uint64 c = c_vector[i].s0;
	if (c == 0) return;
	const uint64 a2k = a2k_vector[i].s0, a2kq = a2k * q;

	for (size_t i = 0; i < factors_loop; ++i)
	{
		if (c % 2 != 0) c = p - c;
		if (c <= 2000000000)
		{
			const uint factor_index = atomic_inc(factor_count);
			factor[factor_index] = (ulong2)(p, c << 32);
		}
		c = mul_mod_const(c, a2k, p, q, a2kq);		// c = a^{(2*i + 1).k}
	}

	c_vector[i].s0 = c;
}

__kernel
void clear_primes(__global uint * restrict const prime_count)
{
	*prime_count = 0;
}
