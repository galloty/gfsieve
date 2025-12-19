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
	uint80 r;
#ifdef PTX_ASM
	const uint32 xs1 = x.s1, ys1 = y.s1;
	uint32 rs1;
	asm volatile ("add.cc.u64 %0, %1, %2;" : "=l" (r.s0) : "l" (x.s0), "l" (y.s0));
	asm volatile ("addc.u32 %0, %1, %2;" : "=r" (rs1) : "r" (xs1), "r" (ys1));
	r.s1 = (uint16)(rs1);
#else
	r.s0 = x.s0 + y.s0; r.s1 = x.s1 + y.s1 + ((r.s0 < x.s0) ? 1 : 0);
#endif
	return r;
}

inline uint80 sub80(const uint80 x, const uint80 y)
{
	uint80 r;
#ifdef PTX_ASM
	const uint32 xs1 = x.s1, ys1 = y.s1;
	uint32 rs1;
	asm volatile ("sub.cc.u64 %0, %1, %2;" : "=l" (r.s0) : "l" (x.s0), "l" (y.s0));
	asm volatile ("subc.u32 %0, %1, %2;" : "=r" (rs1) : "r" (xs1), "r" (ys1));
	r.s1 = (uint16)(rs1);
#else
	 r.s0 = x.s0 - y.s0; r.s1 = x.s1 - y.s1 - ((x.s0 < y.s0) ? 1 : 0);
#endif
	return r;
}

inline uint80 neg80(const uint80 x)
{
	uint80 r;
#ifdef PTX_ASM
	const uint32 xs1 = x.s1;
	uint32 rs1;
	asm volatile ("sub.cc.u64 %0, 0, %1;" : "=l" (r.s0) : "l" (x.s0));
	asm volatile ("subc.u32 %0, 0, %1;" : "=r" (rs1) : "r" (xs1));
	r.s1 = (uint16)(rs1);
#else
	r.s0 = -x.s0; r.s1 = -x.s1 - ((x.s0 != 0) ? 1 : 0);
#endif
	return r;
}

// 0 < s < 64
inline uint80 shl80(const uint80 x, const int s)
{
	const uint16 rs1 = (s < 16) ? (x.s1 << s) : 0;
	uint80 r; r.s0 = x.s0 << s; r.s1 = rs1 | (uint16)(x.s0 >> (64 - s));
	return r;
}

// 0 < s < 64
inline uint80 shr80(const uint80 x, const int s)
{
	const uint16 rs1 = (s < 16) ? (x.s1 >> s) : 0;
	uint80 r; r.s0 = (x.s0 >> s) | ((uint64)(x.s1) << (64 - s)); r.s1 = rs1;
	return r;
}

typedef struct _uint96
{
	uint64 s0;
	uint32 s1;
} uint96;

inline uint96 madd96(const uint96 z, const uint64 x, const uint32 y)
{
	uint96 r;
#ifdef PTX_ASM
	const uint32 xl = lo32(x), xh = hi32(x);
	uint32 c0 = lo32(z.s0), c1 = hi32(z.s0), c2 = z.s1;
	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c0) : "r" (xl), "r" (y), "r" (c0));
	asm volatile ("madc.hi.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (xl), "r" (y), "r" (c1));
	asm volatile ("addc.u32 %0, %1, 0;" : "=r" (c2) : "r" (c2));
	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (xh), "r" (y), "r" (c1));
	asm volatile ("madc.hi.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (xh), "r" (y), "r" (c2));
	r.s0 = upsample(c1, c0); r.s1 = c2;
#else
	r.s0 = z.s0 + x * y; r.s1 = z.s1 + (uint32)(mul_hi(x, (uint64)(y))) + ((r.s0 < z.s0) ? 1 : 0);
#endif
	return r;
}

inline void mul80_wide(const uint80 x, const uint80 y, uint80 * const lo, uint80 * const hi)
{
	lo->s0 = x.s0 * y.s0;
	uint96 t; t.s0 = mul_hi(x.s0, y.s0); t.s1 = x.s1 * (uint32)(y.s1);
	const uint96 r = madd96(madd96(t, x.s0, y.s1), y.s0, x.s1);
	lo->s1 = (uint16)(r.s0); hi->s0 = (r.s0 >> 16) | ((uint64)(r.s1) << 48); hi->s1 = (uint16)(r.s1 >> 16);
}

inline void sqr80_wide(const uint80 x, uint80 * const lo, uint80 * const hi)
{
	lo->s0 = x.s0 * x.s0;
	uint96 t; t.s0 = mul_hi(x.s0, x.s0); t.s1 = x.s1 * (uint32)(x.s1);
	const uint96 r = madd96(t, x.s0, (uint32)(x.s1) << 1);
	lo->s1 = (uint16)(r.s0); hi->s0 = (r.s0 >> 16) | ((uint64)(r.s1) << 48); hi->s1 = (uint16)(r.s1 >> 16);
}

inline uint80 mul80(const uint80 x, const uint80 y)
{
	const uint32 a0 = (uint32)(x.s0), a1 = (uint32)(x.s0 >> 32), a2 = x.s1;
	const uint32 b0 = (uint32)(y.s0), b1 = (uint32)(y.s0 >> 32), b2 = y.s1;
	uint32 c0 = a0 * b0, c1 = mul_hi(a0, b0), c2 = a1 * b1;
#ifdef PTX_ASM
	asm volatile ("mad.lo.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a0), "r" (b2), "r" (c2));
	asm volatile ("mad.lo.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a2), "r" (b0), "r" (c2));
	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a0), "r" (b1), "r" (c1));
	asm volatile ("madc.hi.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a0), "r" (b1), "r" (c2));
	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a1), "r" (b0), "r" (c1));
	asm volatile ("madc.hi.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a1), "r" (b0), "r" (c2));
#else
	c2 += a0 * b2 + a2 * b0;
	const uint64 c12 = upsample(c2, c1) + a0 * (uint64)(b1) + a1 * (uint64)(b0);
	c1 = lo32(c12); c2 = hi32(c12);
#endif
	uint80 r; r.s0 = upsample(c1, c0); r.s1 = (uint16)(c2);
	return r;
}

inline uint80 mul80_hi(const uint80 x, const uint80 y)
{
	uint96 t; t.s0 = mul_hi(x.s0, y.s0); t.s1 = x.s1 * (uint32)(y.s1);
	const uint96 r96 = madd96(madd96(t, x.s0, y.s1), y.s0, x.s1);
	uint80 r; r.s0 = (r96.s0 >> 16) | ((uint64)(r96.s1) << 48); r.s1 = (uint16)(r96.s1 >> 16);
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

inline uint160 sub160(const uint160 x, const uint160 y)
{
	uint160 r;
#ifdef PTX_ASM
	asm volatile ("sub.cc.u64 %0, %1, %2;" : "=l" (r.s0) : "l" (x.s0), "l" (y.s0));
	asm volatile ("subc.cc.u64 %0, %1, %2;" : "=l" (r.s1) : "l" (x.s1), "l" (y.s1));
	asm volatile ("subc.u32 %0, %1, %2;" : "=r" (r.s2) : "r" (x.s2), "r" (y.s2));
#else
	const uint32 c0 = (x.s0 < y.s0) ? 1 : 0, c1 = (x.s1 < y.s1) ? 1 : 0;
	r.s0 = x.s0 - y.s0; r.s1 = x.s1 - y.s1; r.s2 = x.s2 - y.s2;
	const uint32 c2 = (r.s1 < c0) ? 1 : 0;
	r.s1 -= c0; r.s2 -= c1 + c2;
#endif
	return r;
}

// 0 < s < 64
inline uint160 shl160(const uint160 x, const int s)
{
	const uint32 rs2 = (s < 32) ? (x.s2 << s) : 0;
	uint160 r; r.s0 = x.s0 << s; r.s1 = (x.s1 << s) | (x.s0 >> (64 - s)); r.s2 = rs2 | (uint32)(x.s1 >> (64 - s));
	return r;
}

// 0 < s < 64
inline uint160 shr160(const uint160 x, const int s)
{
	const uint32 rs2 = (s < 32) ? (x.s2 >> s) : 0;
	uint160 r; r.s0 = (x.s0 >> s) | (x.s1 << (64 - s)); r.s1 = (x.s1 >> s) | ((uint64)(x.s2) << (64 - s)); r.s2 = rs2;
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

inline uint80 to_int(const uint80 a, const uint80 p, const uint80 q)
{
	const uint80 mp = mul80_hi(mul80(a, q), p);
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
