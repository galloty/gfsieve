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

#include <random>

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
inline bool prp_slow(const __uint128_t p)
{
	mpz_t m, e, r; mpz_inits(m, e, r, nullptr);
	mpz_set_u64(m, uint64_t(p >> 64)); mpz_mul_2exp(m, m, 64); mpz_set_u64(r, uint64_t(p)); mpz_add(m, m, r);
	mpz_sub_ui(e, m, 1); mpz_set_ui(r, 2);
	mpz_powm(r, r, e, m);
	const bool is_prp = (mpz_cmp_ui(r, 1) == 0);
	mpz_clears(m, e, r, nullptr);
	return is_prp;
}
#endif

typedef uint8_t		uint8;
typedef int8_t		int8;
typedef uint16_t	uint16;
typedef uint32_t	uint32;
typedef uint64_t	uint64;

typedef struct _uint80
{
	uint64 s0;
	uint16 s1;
} uint80;

// OpenCL functions
inline uint32 atomic_inc(uint32 * const p) { const uint32 t = *p; (*p)++; return t; }
inline uint32 mul_hi(const uint32 x, const uint32 y) { return uint32((x * uint64(y)) >> 32); }
inline uint64 upsample(const uint32 hi, const uint32 lo) { return ((uint64)(hi) << 32) | lo; }

inline int ilog2_64(const uint64 x) { return 63 - __builtin_clzll((unsigned long long)(x)); }
inline int ilog2_80(const uint80 x) { return (x.s1 != 0) ? (64 + 31 - __builtin_clz((unsigned int)(x.s1))) : ilog2_64(x.s0); }
inline bool bittest_64(const uint64 x, const int b) { return ((x & ((uint64)(1) << b)) != 0); }
inline bool bittest_80(const uint80 x, const int b) { return (b >= 64) ? ((x.s1 & ((uint16)(1) << (b - 64))) != 0) : bittest_64(x.s0, b); }

inline uint32 lo32(const uint64 x) { return (uint32)(x); }
inline uint32 hi32(const uint64 x) { return (uint32)(x >> 32); }

inline uint64 lo32x(const uint64 x) { return (uint64)((uint32)(x)); }
inline uint64 hi32x(const uint64 x) { return x >> 32; }

inline uint32 mul16w(const uint16 x, const uint16 y) { return x * (uint32)(y); }
inline uint64 mul32w(const uint32 x, const uint32 y) { return x * (uint64)(y); }

inline bool eq80(const uint80 x, const uint80 y)	// x is equal to y
{
	return ((x.s1 == y.s1) && (x.s0 == y.s0));
}

inline bool ge80(const uint80 x, const uint80 y)	// x is greater than or equal to y
{
	if (x.s1 > y.s1) return true;
	if (x.s1 < y.s1) return false;
	return (x.s0 >= y.s0);
}

inline bool z80(const uint80 x)	// x is equal to 0
{
	return ((x.s1 == 0) && (x.s0 == 0));
}

inline uint80 add80(const uint80 x, const uint80 y)
{
	uint80 r; r.s0 = x.s0 + y.s0; r.s1 = x.s1 + y.s1 + ((r.s0 < x.s0) ? 1 : 0);
	return r;
}

inline uint80 sub80(const uint80 x, const uint80 y)
{
	uint80 r; r.s0 = x.s0 - y.s0; r.s1 = x.s1 - y.s1 - ((x.s0 < y.s0) ? 1 : 0);
	return r;
}

inline uint80 neg80(const uint80 x)
{
	uint80 r; r.s0 = -x.s0; r.s1 = -x.s1 - ((x.s0 != 0) ? 1 : 0);
	return r;
}

// 0 < s < 64
inline uint80 shl80(const uint80 x, const int s)
{
	const uint16 rs1 = (s < 16) ? (x.s1 << s) : 0;
	uint80 r; r.s1 = rs1 | (uint16)(x.s0 >> (64 - s)); r.s0 = x.s0 << s;
	return r;
}

// 0 < s < 64
inline uint80 shr80(const uint80 x, const int s)
{
	uint80 r; r.s0 = (x.s0 >> s) | ((uint64)(x.s1) << (64 - s)); r.s1 = (s < 16) ? (x.s1 >> s) : 0;
	return r;
}

inline void mul80_wide(const uint80 x, const uint80 y, uint80 & lo, uint80 & hi)
{
	const uint32 a0 = (uint32)(x.s0), a1 = (uint32)(x.s0 >> 32); const uint16 a2 = x.s1;
	const uint32 b0 = (uint32)(y.s0), b1 = (uint32)(y.s0 >> 32); const uint16 b2 = y.s1;

	const uint64 a0b0 = mul32w(a0, b0), a0b1 = mul32w(a0, b1), a1b0 = mul32w(a1, b0), a1b1 = mul32w(a1, b1);	// 64 bits
	const uint64 a0b2 = mul32w(a0, b2), a1b2 = mul32w(a1, b2), a2b0 = mul32w(a2, b0), a2b1 = mul32w(a2, b1);	// 48 bits
	const uint32 a2b2 = mul16w(a2, b2);

	const uint64 c12 = hi32x(a0b0) + lo32x(a0b1) + lo32x(a1b0);
	const uint64 c23 = hi32x(a0b1) + hi32x(a1b0) + lo32x(a1b1) + a0b2 + a2b0 + hi32x(c12);
	const uint64 c34 = hi32x(a1b1) + a1b2 + a2b1 + hi32x(c23);
	const uint32 c0 = lo32(a0b0), c1 = lo32(c12), c2 = lo32(c23), c3 = lo32(c34), c4 = a2b2 + hi32(c34);

	lo.s0 = upsample(c1, c0); lo.s1 = uint16(c2);
	hi.s0 = upsample((c4 << 16) | (c3 >> 16), (c3 << 16) | (c2 >> 16)); hi.s1 = uint16(c4 >> 16);
}

inline void sqr80_wide(const uint80 x, uint80 & lo, uint80 & hi)
{
	const uint32 a0 = (uint32)(x.s0), a1 = (uint32)(x.s0 >> 32); const uint16 a2 = x.s1;

	const uint64 b00 = mul32w(a0, a0), b01 = mul32w(a0, a1), b11 = mul32w(a1, a1);	// 64 bits
	const uint64 b02 = mul32w(a0, a2), b12 = mul32w(a1, a2);	// 48 bits
	const uint32 b22 = mul16w(a2, a2);

	const uint64 c12 = hi32x(b00) + 2 * lo32x(b01);
	const uint64 c23 = 2 * hi32x(b01) + lo32x(b11) + 2 * b02 + hi32x(c12);
	const uint64 c34 = hi32x(b11) + 2 * b12 + hi32x(c23);
	const uint32 c0 = lo32(b00), c1 = lo32(c12), c2 = lo32(c23), c3 = lo32(c34), c4 = b22 + hi32(c34);

	lo.s0 = upsample(c1, c0); lo.s1 = uint16(c2);
	hi.s0 = upsample((c4 << 16) | (c3 >> 16), (c3 << 16) | (c2 >> 16)); hi.s1 = uint16(c4 >> 16);
}

inline uint80 mul80(const uint80 x, const uint80 y)
{
	const uint32 a0 = (uint32)(x.s0), a1 = (uint32)(x.s0 >> 32); const uint16 a2 = x.s1;
	const uint32 b0 = (uint32)(y.s0), b1 = (uint32)(y.s0 >> 32); const uint16 b2 = y.s1;

	const uint64 a0b0 = mul32w(a0, b0), a0b1 = mul32w(a0, b1), a1b0 = mul32w(a1, b0);

	const uint64 c12 = hi32x(a0b0) + a0b1 + a1b0;
	const uint32 c0 = lo32(a0b0), c1 = lo32(c12);
	const uint16 c2 = (uint16)(a1) * (uint16)(b1) + (uint16)(a0) * b2 + a2 * (uint16)(b0) + (uint16)(c12 >> 32);

	uint80 r; r.s0 = upsample(c1, c0); r.s1 = c2;
	return r;
}

inline uint80 mul80_hi(const uint80 x, const uint80 y)
{
	const uint32 a0 = (uint32)(x.s0), a1 = (uint32)(x.s0 >> 32); const uint16 a2 = x.s1;
	const uint32 b0 = (uint32)(y.s0), b1 = (uint32)(y.s0 >> 32); const uint16 b2 = y.s1;

	const uint64 a0b1 = mul32w(a0, b1), a1b0 = mul32w(a1, b0), a1b1 = mul32w(a1, b1);	// 64 bits
	const uint64 a0b2 = mul32w(a0, b2), a1b2 = mul32w(a1, b2), a2b0 = mul32w(a2, b0), a2b1 = mul32w(a2, b1);	// 48 bits
	const uint32 a2b2 = mul16w(a2, b2);

	const uint64 c12 = mul_hi(a0, b0) + lo32x(a0b1) + lo32x(a1b0);
	const uint64 c23 = hi32x(a0b1) + hi32x(a1b0) + lo32x(a1b1) + a0b2 + a2b0 + hi32x(c12);
	const uint64 c34 = hi32x(a1b1) + a1b2 + a2b1 + hi32x(c23);
	const uint32 c2 = lo32(c23), c3 = lo32(c34), c4 = a2b2 + hi32(c34);

	uint80 r; r.s0 = upsample((c4 << 16) | (c3 >> 16), (c3 << 16) | (c2 >> 16)); r.s1 = uint16(c4 >> 16);
	return r;
}

typedef struct _uint160
{
	uint64 s0, s1;
	uint32 s2;
} uint160;

inline bool ge160(const uint160 x, const uint160 y)	// x is greater than or equal to y
{
	if (x.s2 > y.s2) return true;
	if (x.s2 < y.s2) return false;
	if (x.s1 > y.s1) return true;
	if (x.s1 < y.s1) return false;
	return (x.s0 >= y.s0);
}

inline uint160 sub160(const uint160 a, const uint160 b)
{
	const uint32 c0 = (a.s0 < b.s0) ? 1 : 0, c1 = (a.s1 < b.s1) ? 1 : 0;
	uint160 r; r.s0 = a.s0 - b.s0; r.s1 = a.s1 - b.s1; r.s2 = a.s2 - b.s2;
	const uint32 c2 = (r.s1 < c0) ? 1 : 0;
	r.s1 -= c0; r.s2 -= c1 + c2;
	return r;
}

// 0 < s < 64
inline uint160 shl160(const uint160 x, const int s)
{
	const uint32 rs2 = (s < 32) ? (x.s2 << s) : 0;
	uint160 r; r.s2 = rs2 | (uint32)(x.s1 >> (64 - s)); r.s1 = (x.s1 << s) | (x.s0 >> (64 - s)); r.s0 = x.s0 << s;
	return r;
}

// 0 < s < 64
inline uint160 shr160(const uint160 x, const int s)
{
	uint160 r; r.s0 = (x.s0 >> s) | (x.s1 << (64 - s)); r.s1 = (x.s1 >> s) | ((uint64)(x.s2) << (64 - s)); r.s2 = (s < 32) ? (x.s2 >> s) : 0;
	return r;
}

// 2^80 * x / y, x < y, y is not a power of 2
inline uint80 div80(const uint80 x, const uint80 y)
{
	if (z80(x)) return x;

	// 2^b * y <= 2^80 * x < 2^{b+1} * y, 0 <= b < 80
	int b = 80 + ilog2_80(x) - ilog2_80(y) - 1;

	// m = 2^b * y
	const int b_64 = (b < 64) ? 0 : 1, b64 = b % 64;
	uint160 m; m.s0 = (b_64 == 0) ? y.s0 : 0; m.s1 = (b_64 == 0) ? (uint64)(y.s1) : y.s0; m.s2 = (b_64 == 0) ? 0 : (uint32)(y.s1);
	if (b64 > 0) m = shl160(m, b64);

	// z = 2^80 * x
	uint160 z; z.s0 = 0; z.s1 = x.s0 << 16; z.s2 = ((uint32)(x.s1) << 16) | (uint32)(x.s0 >> 48);
	// z < 2*m ?
	uint160 t = shl160(m, 1);
	if (ge160(z, t)) { ++b; m = t; }

	uint80 r; r.s0 = 0; r.s1 = 0;
	for (int j = 0; j <= b; ++j)
	{
		r = shl80(r, 1);
		if (ge160(z, m)) { z = sub160(z, m); r.s0 |= 1; }
		m = shr160(m, 1);
	}

	return r;
}

// 2^80 mod x, x is not a power of 2, 2^63 < x < 2^79
inline uint80 modinv80(const uint80 x)
{
	// 2^b * x < 2^80 < 2^{b+1} * x, 1 <= b <= 16
	const int b = 80 - ilog2_80(x) - 1;

	// m = 2^b * x
	uint80 m = shl80(x, b);
	// r = 2^80 - m
	uint80 r = neg80(m);

	for (int j = 1; j <= b; ++j)
	{
		m = shr80(m, 1);
		if (ge80(r, m)) r = sub80(r, m);
	}

	return r;
}

inline uint80 add_mod(const uint80 a, const uint80 b, const uint80 p)
{
	uint80 r = add80(a, b);
	if (ge80(r, p)) r = sub80(r, p);
	return r;
}

inline uint80 dup_mod(const uint80 a, const uint80 p)
{
	uint80 r = shl80(a, 1);
	if (ge80(r, p)) r = sub80(r, p);
	return r;
}

inline uint80 sub_mod(const uint80 a, const uint80 b, const uint80 p)
{
	uint80 r = sub80(a, b);
	if (!ge80(a, b)) r = add80(r, p);
	return r;
}

// Montgomery modular multiplication
inline uint80 mul_mod(const uint80 a, const uint80 b, const uint80 p, const uint80 q)
{
	uint80 ab_l, ab_h; mul80_wide(a, b, ab_l, ab_h);
	return sub_mod(ab_h, mul80_hi(mul80(ab_l, q), p), p);
}

inline uint80 sqr_mod(const uint80 a, const uint80 p, const uint80 q)
{
	uint80 a2_l, a2_h; sqr80_wide(a, a2_l, a2_h);
	return sub_mod(a2_h, mul80_hi(mul80(a2_l, q), p), p);
}

// Victor Shoup’s modular multiplication
inline uint80 mul_mod_vs(const uint80 a, const uint80 b, const uint80 bp, const uint80 p)
{
	const uint80 abp_h = mul80_hi(a, bp);
	const uint80 r = sub80(mul80(a, b), mul80(abp_h, p));
	return ge80(r, p) ? sub80(r, p) : r;
}

inline uint80 to_int(const uint80 r, const uint80 p, const uint80 q)
{
	const uint80 mp = mul80_hi(mul80(r, q), p);
	return !z80(mp) ? sub80(p, mp) : mp;
}

// p * p_inv = 1 (mod 2^80) (Newton's method)
inline uint80 invert(const uint80 p)
{
	uint80 two; two.s0 = 2; two.s1 = 0;
	uint80 p_inv = sub80(two, p);
	uint80 prev; do { prev = p_inv; p_inv = mul80(p_inv, sub80(two, mul80(p, p_inv))); } while (!eq80(p_inv, prev));
	return p_inv;
}

// a^e mod p, left-to-right algorithm
inline uint80 pow_mod(const uint80 a, const uint64 e, const uint80 p, const uint80 q)
{
	uint80 r = a;
	for (int b = ilog2_64(e) - 1; b >= 0; --b)
	{
		r = sqr_mod(r, p, q);
		if (bittest_64(e, b)) r = mul_mod(r, a, p, q);
	}
	return r;
}

// 2^{(p - 1)/2} ?= +/-1 mod p
inline bool prp(const uint80 p, const uint80 q, const uint80 one)
{
	uint80 e; e.s0 = (p.s0 >> 1) | (uint64(p.s1) << 63); e.s1 = p.s1 >> 1;
	int b = ilog2_80(e) - 1;
	uint80 r = dup_mod(one, p);	// 2 = 1 + 1
	r = dup_mod(r, p);			// 2 * 2 = 2 + 2
	if (bittest_80(e, b)) r = dup_mod(r, p);
	for (--b; b >= 0; --b)
	{
		r = sqr_mod(r, p, q);
		if (bittest_80(e, b)) r = dup_mod(r, p);
	}
	return (eq80(r, one) || eq80(r, sub80(p, one)));
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

static void test(const uint64 i_min, const uint64 i_max, const int n, const int log2_global_worksize, const uint8 * const wheel, const int8 * const kro_vector)
{
	const size_t global_worksize = size_t(1) << log2_global_worksize;
	const size_t factors_loop = size_t(1) << 10;
	const size_t N_2_factors_loop = (size_t(1) << (n - 1)) / factors_loop;
	const size_t prime_size = global_worksize;
	const int g_n = n + 1;

	uint64 * const k_vector = new uint64[prime_size];
	uint80 * const q_vector = new uint80[prime_size];
	uint80 * const ext_vector = new uint80[prime_size];
	uint80 * const c_vector = new uint80[prime_size];
#ifdef CHECK
	uint80 * const c0sqm_vector = new uint80[prime_size];
	uint80 * const qs_vector = new uint80[prime_size];
#endif

	uint32 prime_count = 0;

	for (uint64 i = i_min; i <= i_max; ++i)
	{
		// generate_primes
		for (uint64 id = 0; id < global_worksize; ++id)
		{
			const uint64 j = (i << log2_global_worksize) | id;
			const uint64 k = 15 * (j / 8) + wheel[j % 8];

			uint80 p; p.s0 = (k << g_n) | 1; p.s1 = uint16(k >> (64 - g_n));
			const uint80 q = invert(p);
			// one is the Montgomery form of 1: 2^80 mod p
			const uint80 one = modinv80(p);

			const bool isprp = prp(p, q, one);
#ifdef CHECK
			if (isprp != prp_slow(p.s0 | (__uint128_t(p.s1) << 64)))
			{
				std::cout << "p = " << p.s0 << "," << p.s1 << std::endl;
				throw std::runtime_error("Error: prp.");
			}
#endif
			if (isprp)
			{
				const uint32 prime_index = atomic_inc(&prime_count);
				k_vector[prime_index] = k;
				q_vector[prime_index] = q;
				ext_vector[prime_index] = one;
			}
		}

		std::cout << i << ": " << prime_count << " primes" << std::endl;

		// init_factors
		for (uint64 id = 0; id < global_worksize; ++id)
		{
			if (id >= prime_count) break;

			const uint64 k = k_vector[id];
			const uint80 q = q_vector[id], one = ext_vector[id];

			uint80 p; p.s0 = (k << g_n) | 1; p.s1 = uint16(k >> (64 - g_n));
			const uint80 two = dup_mod(one, p);

			// p = 1 (mod 4). If a is odd then (a/p) = (p/a) = ({p mod a}/a)

			uint32 a = 3; uint80 am = add_mod(two, one, p);
			const uint32 pmod3 = ((uint32(k % 3) << g_n) | 1) % 3;
			if (pmod3 != 2)
			{
				a += 2; am = add_mod(am, two, p);
				const uint32 pmod5 = ((uint32(k % 5) << g_n) | 1) % 5;
				if (kro_vector[256 * ((5 - 3) / 2) + pmod5] >= 0)
				{
					a += 2; am = add_mod(am, two, p);
					const uint32 pmod7 = ((uint32(k % 7) << g_n) | 1) % 7;
					if (kro_vector[256 * ((7 - 3) / 2) + pmod7] >= 0)
					{
						a += 4; am = add_mod(am, dup_mod(two, p), p);
						while (a < 256)
						{
							const uint32 pmoda = uint32((((k % a) << g_n) | 1) % a);
							if (kro_vector[256 * ((a - 3) / 2) + pmoda] < 0) break;
							a += 2; am = add_mod(am, two, p);
						}
						if (a >= 256) { a = 0; am.s0 = 0; am.s1 = 0; }
					}
				}
			}

			// a^{(p - 1)/2} = -1 <=> a^{k*2^n} = -1. (a^k)^{2*i + 1} are the roots of b^{2^n} + 1 = 0 (mod p)
			const uint80 cm = pow_mod(am, k, p, q);
#ifdef CHECK
			// p - 1 = k * 2^n
			uint80 t = cm; for (int i = 0; i < n; ++i) t = sqr_mod(t, p, q);
			if (!eq80(t, sub80(p, one))) throw std::runtime_error("Error: a.");
#endif
			const uint80 c = to_int(cm, p, q), c2 = mul_mod(c, cm, p, q);
			c_vector[id] = c;
			ext_vector[id] = c2;
			q_vector[id] = div80(c2, p);
#ifdef CHECK
			const uint80 cm2 = sqr_mod(cm, p, q);
			c0sqm_vector[id] = cm2;
			qs_vector[id] = q;
#endif
		}

		// check_factors
		for (size_t j = 0; j < N_2_factors_loop; ++j)
		{
			std::cout << j << " / " << N_2_factors_loop << std::endl;

			for (uint64 id = 0; id < global_worksize; ++id)
			{
				if (id >= prime_count) break;

				uint80 c = c_vector[id];
				if (z80(c)) continue;

				const uint64 k = k_vector[id];
				uint80 p; p.s0 = (k << g_n) | 1; p.s1 = uint32(k >> (64 - g_n));

				const uint80 c0sq = ext_vector[id], c0sqp = q_vector[id];
#ifdef CHECK
				const uint80 c0sqm = c0sqm_vector[id], q = qs_vector[id];
#endif
				for (size_t l = 0; l < factors_loop; ++l)
				{
					const uint80 b = (c.s0 % 2 == 0) ? c : sub80(p, c);
					const bool found = ((b.s1 == 0) && (b.s0 <= 2000000000));

					if (found) check(k, uint32(b.s0), n);
#ifdef CHECK
					const uint80 cp = c;
#endif
					c = mul_mod_vs(c, c0sq, c0sqp, p);	// c = a^{(2*i + 1).k}
#ifdef CHECK
					if (!eq80(c, mul_mod(cp, c0sqm, p, q))) throw std::runtime_error("Error: factors loop.");
#endif
				}

				c_vector[id] = c;
			}
		}

		prime_count = 0;
	}

	delete[] k_vector;
	delete[] q_vector;
	delete[] ext_vector;
	delete[] c_vector;
#ifdef CHECK
	delete[] c0sqm_vector;
	delete[] qs_vector;
#endif
}

int main()
{
	try
	{
		// Kronecker symbols (i/j) for odd j <= 255.
		int8 kro_vector[128 * 256];
		mpz_t zj; mpz_init(zj);
		for (uint32_t j = 3; j < 256; j += 2)
		{
			mpz_set_ui(zj, j);
			for (uint32_t i = 0; i < j; ++i)
			{
				kro_vector[256 * ((j - 3) / 2) + i] = int8(mpz_ui_kronecker(i, zj));
			}
		}
		mpz_clear(zj);

		const int n = 23;

		// gcd(p mod 15, 15) = 1: 8 solutions
		uint8 wheel[8];
		size_t i = 0;
		for (uint64 k = 0; k < 15; ++k)
		{
			const uint64 p = (k << (n + 1)) | 1;
			if ((p % 3 == 0) || (p % 5 == 0)) continue;
			wheel[i] = uint8(k);
			++i;
		}
		if (i != 8) throw std::runtime_error("Error: wheel count.");

		const int log2_global_worksize = 22;
		const double unit = 1e15;

		const uint32 p_min = 20000, p_max = p_min + 1;	// 18446 <= p <= 600000000

		const double f = unit * 8.0 / 15 / std::pow(2.0, double(n + 1 + log2_global_worksize));
		const uint64 i_min = uint64(floor(p_min * f)), i_max = uint64(ceil(p_max * f));

		mpz_t zp_min, zp_max; mpz_inits(zp_min, zp_max, nullptr);

		const uint64 j_min = uint64(i_min) << log2_global_worksize, j_max = uint64(i_max) << log2_global_worksize;
		const uint64 k_min = 15 * (j_min / 8) + wheel[j_min % 8], k_max = 15 * (j_max / 8) + wheel[j_max % 8];

		mpz_set_u64(zp_min, k_min); mpz_mul_2exp(zp_min, zp_min, n + 1); mpz_add_ui(zp_min, zp_min, 1);
		mpz_set_u64(zp_max, k_max); mpz_mul_2exp(zp_max, zp_max, n + 1); mpz_add_ui(zp_max, zp_max, 1);

		char p_min_str[32], p_max_str[32]; mpz_get_str(p_min_str, 10, zp_min); mpz_get_str(p_max_str, 10, zp_max);
		std::cout << "Testing n = " << n << ", p in [" << p_min_str << ", " << p_max_str << "[." << std::endl;

		mpz_clears(zp_min, zp_max, nullptr);

		std::cout << "For i = " << i_min << " to " << i_max << std::endl;

		test(i_min, i_max, n, log2_global_worksize, wheel, kro_vector);
	}
	catch (const std::runtime_error & e)
	{
		std::cout << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}