/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#ifdef __NV_CL_C_VERSION
	#define PTX_ASM	1
#endif

#define sz_t		uint
#ifndef uint_8
#define uint_8		uchar
#endif
#define int_8		char
#define uint16		ushort
#define uint32		uint
#define uint64		ulong
#define uint64_2	ulong2

#if !defined(G_N)
#define G_N			16
#define LN_BLKSZ	22
#define FBLK		1024

__constant uint_8 wheel[8] = { 0, 1, 3, 6, 7, 10, 12, 13 };
#endif

typedef struct _uint80
{
	uint64 s0;
	uint16 s1;
} uint80;

inline int ilog2_64(const uint64 x) { return 63 - clz(x); }
inline int ilog2_80(const uint80 x) { return (x.s1 != 0) ? (64 + 31 - clz((uint32)x.s1)) : ilog2_64(x.s0); }
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

inline void mul80_wide(const uint80 x, const uint80 y, uint80 * const lo, uint80 * const hi)
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

	lo->s0 = upsample(c1, c0); lo->s1 = (uint16)(c2);
	hi->s0 = upsample((c4 << 16) | (c3 >> 16), (c3 << 16) | (c2 >> 16)); hi->s1 = (uint16)(c4 >> 16);
}

inline void sqr80_wide(const uint80 x, uint80 * const lo, uint80 * const hi)
{
	const uint32 a0 = (uint32)(x.s0), a1 = (uint32)(x.s0 >> 32); const uint16 a2 = x.s1;

	const uint64 b00 = mul32w(a0, a0), b01 = mul32w(a0, a1), b11 = mul32w(a1, a1);	// 64 bits
	const uint64 b02 = mul32w(a0, a2), b12 = mul32w(a1, a2);	// 48 bits
	const uint32 b22 = mul16w(a2, a2);

	const uint64 c12 = hi32x(b00) + 2 * lo32x(b01);
	const uint64 c23 = 2 * hi32x(b01) + lo32x(b11) + 2 * b02 + hi32x(c12);
	const uint64 c34 = hi32x(b11) + 2 * b12 + hi32x(c23);
	const uint32 c0 = lo32(b00), c1 = lo32(c12), c2 = lo32(c23), c3 = lo32(c34), c4 = b22 + hi32(c34);

	lo->s0 = upsample(c1, c0); lo->s1 = (uint16)(c2);
	hi->s0 = upsample((c4 << 16) | (c3 >> 16), (c3 << 16) | (c2 >> 16)); hi->s1 = (uint16)(c4 >> 16);
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

	uint80 r; r.s0 = upsample((c4 << 16) | (c3 >> 16), (c3 << 16) | (c2 >> 16)); r.s1 = (uint16)(c4 >> 16);
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
	uint80 ab_l, ab_h; mul80_wide(a, b, &ab_l, &ab_h);
	return sub_mod(ab_h, mul80_hi(mul80(ab_l, q), p), p);
}

inline uint80 sqr_mod(const uint80 a, const uint80 p, const uint80 q)
{
	uint80 a2_l, a2_h; sqr80_wide(a, &a2_l, &a2_h);
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
	uint80 e; e.s0 = (p.s0 >> 1) | ((uint64)(p.s1) << 63); e.s1 = p.s1 >> 1;
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

/*

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
}*/

__kernel
void generate_primes(__global sz_t * restrict const prime_count,
	__global uint64 * restrict const k_vector, __global uint64_2 * restrict const q_vector,
	__global uint64_2 * restrict const ext_vector, const uint64 i)
{
	const sz_t id = (sz_t)get_global_id(0);

	const uint64 j = (i << LN_BLKSZ) | id;
	const uint64 k = 15 * (j / 8) + wheel[j % 8];

	uint80 p; p.s0 = (k << G_N) | 1; p.s1 = (uint16)(k >> (64 - G_N));
	const uint80 q = invert(p);
	// one is the Montgomery form of 1: 2^80 mod p
	const uint80 one = modinv80(p);

	if (prp(p, q, one))
	{
		const uint prime_index = atomic_inc(prime_count);
		k_vector[prime_index] = k;
		q_vector[prime_index] = (uint64_2)(q.s0, q.s1);
		ext_vector[prime_index] = (uint64_2)(one.s0, one.s1);
	}
}

__kernel
void init_factors(__global const sz_t * restrict const prime_count,
	__global const uint64 * restrict const k_vector, __global uint64_2 * restrict const q_vector,
	__global uint64_2 * restrict const ext_vector, __global const int_8 * restrict const kro_vector,
	__global uint64_2 * restrict const c_vector)
{
	const sz_t id = (sz_t)get_global_id(0);
	if (id >= *prime_count) return;

	const uint64 k = k_vector[id];
	const uint64_2 ql = q_vector[id], onel = ext_vector[id];
	uint80 q; q.s0 = ql.s0; q.s1 = (uint16)(ql.s1);
	uint80 one; one.s0 = onel.s0; one.s1 = (uint16)(onel.s1);

	uint80 p; p.s0 = (k << G_N) | 1; p.s1 = (uint16)(k >> (64 - G_N));
	const uint80 two = dup_mod(one, p);

	// p = 1 (mod 4). If a is odd then (a/p) = (p/a) = ({p mod a}/a)

	uint32 a = 3; uint80 am = add_mod(two, one, p);
	const uint32 pmod3 = (((uint32)(k % 3) << G_N) | 1) % 3;
	if (pmod3 != 2)
	{
		a += 2; am = add_mod(am, two, p);
		const uint32 pmod5 = (((uint32)(k % 5) << G_N) | 1) % 5;
		if (kro_vector[256 * ((5 - 3) / 2) + pmod5] >= 0)
		{
			a += 2; am = add_mod(am, two, p);
			const uint32 pmod7 = (((uint32)(k % 7) << G_N) | 1) % 7;
			if (kro_vector[256 * ((7 - 3) / 2) + pmod7] >= 0)
			{
				a += 4; am = add_mod(am, dup_mod(two, p), p);
				while (a < 256)
				{
					const uint32 pmoda = (uint32)((((k % a) << G_N) | 1) % a);
					if (kro_vector[256 * ((a - 3) / 2) + pmoda] < 0) break;
					a += 2; am = add_mod(am, two, p);
				}
				if (a >= 256)
				{
					c_vector[id] = (uint64_2)(0, 0);
					return;
				}
			}
		}
	}

	// a^{(p - 1)/2} = -1 <=> a^{k*2^n} = -1. (a^k)^{2*i + 1} are the roots of b^{2^n} + 1 = 0 (mod p)
	const uint80 cm = pow_mod(am, k, p, q);
	const uint80 c = to_int(cm, p, q), c2 = mul_mod(c, cm, p, q);
	c_vector[id] = (uint64_2)(c.s0, c.s1);
	ext_vector[id] = (uint64_2)(c2.s0, c2.s1);
	const uint80 c2p = div80(c2, p);
	q_vector[id] = (uint64_2)(c2p.s0, c2p.s1);
}

__kernel
void check_factors(__global const sz_t * restrict const prime_count,
	__global const uint64 * restrict const k_vector, __global const uint64_2 * restrict const q_vector,
	__global uint64_2 * restrict const c_vector, __global const uint64_2 * restrict const ext_vector,
	__global sz_t * restrict const factor_count, __global uint64_2 * restrict const factor_vector)
{
	const sz_t id = (sz_t)get_global_id(0);
	if (id >= *prime_count) return;

	const uint64_2 cl = c_vector[id];
	uint80 c; c.s0 = cl.s0; c.s1 = (uint16)(cl.s1);
	if (z80(c)) return;

	const uint64 k = k_vector[id];
	uint80 p; p.s0 = (k << G_N) | 1; p.s1 = (uint16)(k >> (64 - G_N));

	const uint64_2 c0sql = ext_vector[id], c0sqpl = q_vector[id];
	uint80 c0sq; c0sq.s0 = c0sql.s0; c0sq.s1 = (uint16)(c0sql.s1);
	uint80 c0sqp; c0sqp.s0 = c0sqpl.s0; c0sqp.s1 = (uint16)(c0sqpl.s1);

	for (size_t l = 0; l < FBLK; ++l)
	{
		const uint80 b = (c.s0 % 2 == 0) ? c : sub80(p, c);

		if ((b.s1 == 0) && (b.s0 <= 2000000000))
		{
			const sz_t factor_index = atomic_inc(factor_count);
			factor_vector[factor_index] = (uint64_2)(k, b.s0);
		}

		c = mul_mod_vs(c, c0sq, c0sqp, p);	// c = a^{(2*i + 1).k}
	}

	c_vector[id] = (uint64_2)(c.s0, c.s1);
}

__kernel
void clear(__global sz_t * restrict const count)
{
	*count = 0;
}
