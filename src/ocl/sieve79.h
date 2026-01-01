/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>

static const char * const src_ocl_sieve79 = \
"/*\n" \
"Copyright 2020, Yves Gallot\n" \
"\n" \
"gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.\n" \
"Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.\n" \
"*/\n" \
"\n" \
"#ifdef __NV_CL_C_VERSION\n" \
"	#define PTX_ASM	1\n" \
"#endif\n" \
"\n" \
"#if !defined(G_N)\n" \
"#define G_N			16\n" \
"#define LN_BLKSZ	22\n" \
"#define FBLK		1024\n" \
"typedef	uchar		uint_8;\n" \
"\n" \
"__constant uint_8 wheel[8] = { 0, 1, 3, 6, 7, 10, 12, 13 };\n" \
"#endif\n" \
"\n" \
"typedef uint	sz_t;\n" \
"typedef char	int_8;\n" \
"typedef ushort	uint_16;\n" \
"typedef uint	uint_32;\n" \
"typedef ulong	uint_64;\n" \
"typedef ulong2	uint_64_2;\n" \
"\n" \
"typedef struct _uint_80\n" \
"{\n" \
"	uint_64 lo;\n" \
"	uint_16 hi;\n" \
"} uint_80;\n" \
"\n" \
"inline int ilog2_64(const uint_64 x) { return 63 - clz(x); }\n" \
"inline int ilog2_80(const uint_80 x) { return (x.hi != 0) ? (64 + 31 - clz((uint_32)x.hi)) : ilog2_64(x.lo); }\n" \
"inline bool bittest_64(const uint_64 x, const int b) { return ((x & ((uint_64)(1) << b)) != 0); }\n" \
"inline bool bittest_80(const uint_80 x, const int b) { return (b >= 64) ? ((x.hi & ((uint_16)(1) << (b - 64))) != 0) : bittest_64(x.lo, b); }\n" \
"\n" \
"inline uint_32 lo32(const uint_64 x) { return (uint_32)(x); }\n" \
"inline uint_32 hi32(const uint_64 x) { return (uint_32)(x >> 32); }\n" \
"\n" \
"inline bool eq80(const uint_80 x, const uint_80 y)	// x is equal to y\n" \
"{\n" \
"	return ((x.hi == y.hi) && (x.lo == y.lo));\n" \
"}\n" \
"\n" \
"inline bool ge80(const uint_80 x, const uint_80 y)	// x is greater than or equal to y\n" \
"{\n" \
"	if (x.hi > y.hi) return true;\n" \
"	if (x.hi < y.hi) return false;\n" \
"	return (x.lo >= y.lo);\n" \
"}\n" \
"\n" \
"inline bool z80(const uint_80 x)	// x is equal to 0\n" \
"{\n" \
"	return ((x.hi == 0) && (x.lo == 0));\n" \
"}\n" \
"\n" \
"inline uint_80 add80(const uint_80 x, const uint_80 y)\n" \
"{\n" \
"	uint_80 r;\n" \
"#ifdef PTX_ASM\n" \
"	const uint32 xhi = x.hi, yhi = y.hi;\n" \
"	uint32 rhi;\n" \
"	asm volatile (\"add.cc.u64 %0, %1, %2;\" : \"=l\" (r.lo) : \"l\" (x.lo), \"l\" (y.lo));\n" \
"	asm volatile (\"addc.u32 %0, %1, %2;\" : \"=r\" (rhi) : \"r\" (xhi), \"r\" (yhi));\n" \
"	r.hi = (uint_16)(rhi);\n" \
"#else\n" \
"	r.lo = x.lo + y.lo; r.hi = x.hi + y.hi + ((r.lo < x.lo) ? 1 : 0);\n" \
"#endif\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint_80 sub80(const uint_80 x, const uint_80 y)\n" \
"{\n" \
"	uint_80 r;\n" \
"#ifdef PTX_ASM\n" \
"	const uint32 xhi = x.hi, yhi = y.hi;\n" \
"	uint32 rhi;\n" \
"	asm volatile (\"sub.cc.u64 %0, %1, %2;\" : \"=l\" (r.lo) : \"l\" (x.lo), \"l\" (y.lo));\n" \
"	asm volatile (\"subc.u32 %0, %1, %2;\" : \"=r\" (rhi) : \"r\" (xhi), \"r\" (yhi));\n" \
"	r.hi = (uint_16)(rhi);\n" \
"#else\n" \
"	r.lo = x.lo - y.lo; r.hi = x.hi - y.hi - ((x.lo < y.lo) ? 1 : 0);\n" \
"#endif\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint_80 neg80(const uint_80 x)\n" \
"{\n" \
"	uint_80 r;\n" \
"#ifdef PTX_ASM\n" \
"	const uint32 xhi = x.hi;\n" \
"	uint32 rhi;\n" \
"	asm volatile (\"sub.cc.u64 %0, 0, %1;\" : \"=l\" (r.lo) : \"l\" (x.lo));\n" \
"	asm volatile (\"subc.u32 %0, 0, %1;\" : \"=r\" (rhi) : \"r\" (xhi));\n" \
"	r.hi = (uint_16)(rhi);\n" \
"#else\n" \
"	r.lo = -x.lo; r.hi = -x.hi - ((x.lo != 0) ? 1 : 0);\n" \
"#endif\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 0 < s < 64\n" \
"inline uint_80 shl80(const uint_80 x, const int s)\n" \
"{\n" \
"	const uint_16 rhi = (s < 16) ? (x.hi << s) : 0;\n" \
"	uint_80 r; r.lo = x.lo << s; r.hi = rhi | (uint_16)(x.lo >> (64 - s));\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 0 < s < 64\n" \
"inline uint_80 shr80(const uint_80 x, const int s)\n" \
"{\n" \
"	const uint_16 rhi = (s < 16) ? (x.hi >> s) : 0;\n" \
"	uint_80 r; r.lo = (x.lo >> s) | ((uint_64)(x.hi) << (64 - s)); r.hi = rhi;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"typedef struct _uint_96\n" \
"{\n" \
"	uint_64 lo;\n" \
"	uint_32 hi;\n" \
"} uint_96;\n" \
"\n" \
"inline uint_96 madd96(const uint_96 z, const uint_64 x, const uint_32 y)\n" \
"{\n" \
"	uint_96 r;\n" \
"#ifdef PTX_ASM\n" \
"	const uint32 xl = lo32(x), xh = hi32(x);\n" \
"	uint32 c0 = lo32(z.s0), c1 = hi32(z.s0), c2 = z.s1;\n" \
"	asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c0) : \"r\" (xl), \"r\" (y), \"r\" (c0));\n" \
"	asm volatile (\"madc.hi.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c1) : \"r\" (xl), \"r\" (y), \"r\" (c1));\n" \
"	asm volatile (\"addc.u32 %0, %1, 0;\" : \"=r\" (c2) : \"r\" (c2));\n" \
"	asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c1) : \"r\" (xh), \"r\" (y), \"r\" (c1));\n" \
"	asm volatile (\"madc.hi.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (xh), \"r\" (y), \"r\" (c2));\n" \
"	r.s0 = upsample(c1, c0); r.s1 = c2;\n" \
"#else\n" \
"	r.lo = z.lo + x * y; r.hi = z.hi + (uint_32)(mul_hi(x, (uint_64)(y))) + ((r.lo < z.lo) ? 1 : 0);\n" \
"#endif\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline void mul80_wide(const uint_80 x, const uint_80 y, uint_80 * const lo, uint_80 * const hi)\n" \
"{\n" \
"	lo->lo = x.lo * y.lo;\n" \
"	uint_96 t; t.lo = mul_hi(x.lo, y.lo); t.hi = x.hi * (uint_32)(y.hi);\n" \
"	const uint_96 r = madd96(madd96(t, x.lo, y.hi), y.lo, x.hi);\n" \
"	lo->hi = (uint_16)(r.lo); hi->lo = (r.lo >> 16) | ((uint_64)(r.hi) << 48); hi->hi = (uint_16)(r.hi >> 16);\n" \
"}\n" \
"\n" \
"inline void sqr80_wide(const uint_80 x, uint_80 * const lo, uint_80 * const hi)\n" \
"{\n" \
"	lo->lo = x.lo * x.lo;\n" \
"	uint_96 t; t.lo = mul_hi(x.lo, x.lo); t.hi = x.hi * (uint_32)(x.hi);\n" \
"	const uint_96 r = madd96(t, x.lo, (uint_32)(x.hi) << 1);\n" \
"	lo->hi = (uint_16)(r.lo); hi->lo = (r.lo >> 16) | ((uint_64)(r.hi) << 48); hi->hi = (uint_16)(r.hi >> 16);\n" \
"}\n" \
"\n" \
"inline uint_80 mul80(const uint_80 x, const uint_80 y)\n" \
"{\n" \
"	const uint_32 a0 = lo32(x.lo), a1 = hi32(x.lo), a2 = x.hi;\n" \
"	const uint_32 b0 = lo32(y.lo), b1 = hi32(y.lo), b2 = y.hi;\n" \
"	uint_32 c0 = a0 * b0, c1 = mul_hi(a0, b0), c2 = a1 * b1;\n" \
"#ifdef PTX_ASM\n" \
"	asm volatile (\"mad.lo.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a0), \"r\" (b2), \"r\" (c2));\n" \
"	asm volatile (\"mad.lo.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a2), \"r\" (b0), \"r\" (c2));\n" \
"	asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c1) : \"r\" (a0), \"r\" (b1), \"r\" (c1));\n" \
"	asm volatile (\"madc.hi.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a0), \"r\" (b1), \"r\" (c2));\n" \
"	asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c1) : \"r\" (a1), \"r\" (b0), \"r\" (c1));\n" \
"	asm volatile (\"madc.hi.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a1), \"r\" (b0), \"r\" (c2));\n" \
"#else\n" \
"	c2 += a0 * b2 + a2 * b0;\n" \
"	const uint_64 c12 = upsample(c2, c1) + a0 * (uint_64)(b1) + a1 * (uint_64)(b0);\n" \
"	c1 = lo32(c12); c2 = hi32(c12);\n" \
"#endif\n" \
"	uint_80 r; r.lo = upsample(c1, c0); r.hi = (uint_16)(c2);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint_80 mul80_hi(const uint_80 x, const uint_80 y)\n" \
"{\n" \
"	uint_96 t; t.lo = mul_hi(x.lo, y.lo); t.hi = x.hi * (uint_32)(y.hi);\n" \
"	const uint_96 r96 = madd96(madd96(t, x.lo, y.hi), y.lo, x.hi);\n" \
"	uint_80 r; r.lo = (r96.lo >> 16) | ((uint_64)(r96.hi) << 48); r.hi = (uint_16)(r96.hi >> 16);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"typedef struct _uint_160\n" \
"{\n" \
"	uint_64 a0, a1;\n" \
"	uint_32 a2;\n" \
"} uint_160;\n" \
"\n" \
"inline bool ge160(const uint_160 x, const uint_160 y)	// x is greater than or equal to y\n" \
"{\n" \
"	if (x.a2 > y.a2) return true;\n" \
"	if (x.a2 < y.a2) return false;\n" \
"	if (x.a1 > y.a1) return true;\n" \
"	if (x.a1 < y.a1) return false;\n" \
"	return (x.a0 >= y.a0);\n" \
"}\n" \
"\n" \
"inline uint_160 sub160(const uint_160 x, const uint_160 y)\n" \
"{\n" \
"	uint_160 r;\n" \
"#ifdef PTX_ASM\n" \
"	asm volatile (\"sub.cc.u64 %0, %1, %2;\" : \"=l\" (r.a0) : \"l\" (x.a0), \"l\" (y.a0));\n" \
"	asm volatile (\"subc.cc.u64 %0, %1, %2;\" : \"=l\" (r.a1) : \"l\" (x.a1), \"l\" (y.a1));\n" \
"	asm volatile (\"subc.u32 %0, %1, %2;\" : \"=r\" (r.a2) : \"r\" (x.a2), \"r\" (y.a2));\n" \
"#else\n" \
"	const uint_32 c0 = (x.a0 < y.a0) ? 1 : 0, c1 = (x.a1 < y.a1) ? 1 : 0;\n" \
"	r.a0 = x.a0 - y.a0; r.a1 = x.a1 - y.a1; r.a2 = x.a2 - y.a2;\n" \
"	const uint_32 c2 = (r.a1 < c0) ? 1 : 0;\n" \
"	r.a1 -= c0; r.a2 -= c1 + c2;\n" \
"#endif\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 0 < s < 64\n" \
"inline uint_160 shl160(const uint_160 x, const int s)\n" \
"{\n" \
"	const uint_32 ra2 = (s < 32) ? (x.a2 << s) : 0;\n" \
"	uint_160 r; r.a0 = x.a0 << s; r.a1 = (x.a1 << s) | (x.a0 >> (64 - s)); r.a2 = ra2 | (uint_32)(x.a1 >> (64 - s));\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 0 < s < 64\n" \
"inline uint_160 shr160(const uint_160 x, const int s)\n" \
"{\n" \
"	const uint_32 ra2 = (s < 32) ? (x.a2 >> s) : 0;\n" \
"	uint_160 r; r.a0 = (x.a0 >> s) | (x.a1 << (64 - s)); r.a1 = (x.a1 >> s) | ((uint_64)(x.a2) << (64 - s)); r.a2 = ra2;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 2^80 * x / y, x < y, y is not a power of 2\n" \
"inline uint_80 div80(const uint_80 x, const uint_80 y)\n" \
"{\n" \
"	if (z80(x)) return x;\n" \
"\n" \
"	// 2^b * y <= 2^80 * x < 2^{b+1} * y, 0 <= b < 80\n" \
"	int b = 80 + ilog2_80(x) - ilog2_80(y) - 1;\n" \
"\n" \
"	// m = 2^b * y\n" \
"	const int b_64 = (b < 64) ? 0 : 1, b64 = b % 64;\n" \
"	uint_160 m; m.a0 = (b_64 == 0) ? y.lo : 0; m.a1 = (b_64 == 0) ? (uint_64)(y.hi) : y.lo; m.a2 = (b_64 == 0) ? 0 : (uint_32)(y.hi);\n" \
"	if (b64 > 0) m = shl160(m, b64);\n" \
"\n" \
"	// z = 2^80 * x\n" \
"	uint_160 z; z.a0 = 0; z.a1 = x.lo << 16; z.a2 = ((uint_32)(x.hi) << 16) | (uint_32)(x.lo >> 48);\n" \
"	// z < 2*m ?\n" \
"	uint_160 t = shl160(m, 1);\n" \
"	if (ge160(z, t)) { ++b; m = t; }\n" \
"\n" \
"	uint_80 r; r.lo = 0; r.hi = 0;\n" \
"	for (int j = 0; j <= b; ++j)\n" \
"	{\n" \
"		r = shl80(r, 1);\n" \
"		if (ge160(z, m)) { z = sub160(z, m); r.lo |= 1; }\n" \
"		m = shr160(m, 1);\n" \
"	}\n" \
"\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 2^80 mod x, x is not a power of 2, 2^63 < x < 2^79\n" \
"inline uint_80 modinv80(const uint_80 x)\n" \
"{\n" \
"	// 2^b * x < 2^80 < 2^{b+1} * x, 1 <= b <= 16\n" \
"	const int b = 80 - ilog2_80(x) - 1;\n" \
"\n" \
"	// m = 2^b * x\n" \
"	uint_80 m = shl80(x, b);\n" \
"	// r = 2^80 - m\n" \
"	uint_80 r = neg80(m);\n" \
"\n" \
"	for (int j = 1; j <= b; ++j)\n" \
"	{\n" \
"		m = shr80(m, 1);\n" \
"		if (ge80(r, m)) r = sub80(r, m);\n" \
"	}\n" \
"\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint_80 add_mod(const uint_80 a, const uint_80 b, const uint_80 p)\n" \
"{\n" \
"	uint_80 r = add80(a, b);\n" \
"	if (ge80(r, p)) r = sub80(r, p);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint_80 dup_mod(const uint_80 a, const uint_80 p)\n" \
"{\n" \
"	uint_80 r = shl80(a, 1);\n" \
"	if (ge80(r, p)) r = sub80(r, p);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint_80 sub_mod(const uint_80 a, const uint_80 b, const uint_80 p)\n" \
"{\n" \
"	uint_80 r = sub80(a, b);\n" \
"	if (!ge80(a, b)) r = add80(r, p);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// Montgomery modular multiplication\n" \
"inline uint_80 mul_mod(const uint_80 a, const uint_80 b, const uint_80 p, const uint_80 q)\n" \
"{\n" \
"	uint_80 ab_l, ab_h; mul80_wide(a, b, &ab_l, &ab_h);\n" \
"	return sub_mod(ab_h, mul80_hi(mul80(ab_l, q), p), p);\n" \
"}\n" \
"\n" \
"inline uint_80 sqr_mod(const uint_80 a, const uint_80 p, const uint_80 q)\n" \
"{\n" \
"	uint_80 a2_l, a2_h; sqr80_wide(a, &a2_l, &a2_h);\n" \
"	return sub_mod(a2_h, mul80_hi(mul80(a2_l, q), p), p);\n" \
"}\n" \
"\n" \
"// Victor Shoup’s modular multiplication\n" \
"inline uint_80 mul_mod_vs(const uint_80 a, const uint_80 b, const uint_80 bp, const uint_80 p)\n" \
"{\n" \
"	const uint_80 abp_h = mul80_hi(a, bp);\n" \
"	const uint_80 r = sub80(mul80(a, b), mul80(abp_h, p));\n" \
"	return ge80(r, p) ? sub80(r, p) : r;\n" \
"}\n" \
"\n" \
"inline uint_80 to_int(const uint_80 a, const uint_80 p, const uint_80 q)\n" \
"{\n" \
"	const uint_80 mr = mul80_hi(mul80(a, q), p);\n" \
"	return !z80(mr) ? sub80(p, mr) : mr;\n" \
"}\n" \
"\n" \
"// p * p_inv = 1 (mod 2^80) (Newton's method)\n" \
"inline uint_80 invert(const uint_80 p)\n" \
"{\n" \
"	uint_80 two; two.lo = 2; two.hi = 0;\n" \
"	uint_80 p_inv = sub80(two, p);\n" \
"	uint_80 prev; do { prev = p_inv; p_inv = mul80(p_inv, sub80(two, mul80(p, p_inv))); } while (!eq80(p_inv, prev));\n" \
"	return p_inv;\n" \
"}\n" \
"\n" \
"// a^e mod p, left-to-right algorithm\n" \
"inline uint_80 pow_mod(const uint_80 a, const uint_64 e, const uint_80 p, const uint_80 q)\n" \
"{\n" \
"	uint_80 r = a;\n" \
"	for (int b = ilog2_64(e) - 1; b >= 0; --b)\n" \
"	{\n" \
"		r = sqr_mod(r, p, q);\n" \
"		if (bittest_64(e, b)) r = mul_mod(r, a, p, q);\n" \
"	}\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 2^{(p - 1)/2} ?= +/-1 mod p\n" \
"inline bool prp(const uint_80 p, const uint_80 q, const uint_80 one)\n" \
"{\n" \
"	uint_80 e; e.lo = (p.lo >> 1) | ((uint_64)(p.hi) << 63); e.hi = p.hi >> 1;\n" \
"	int b = ilog2_80(e) - 1;\n" \
"	uint_80 r = dup_mod(one, p);	// 2 = 1 + 1\n" \
"	r = dup_mod(r, p);			// 2 * 2 = 2 + 2\n" \
"	if (bittest_80(e, b)) r = dup_mod(r, p);\n" \
"	for (--b; b >= 0; --b)\n" \
"	{\n" \
"		r = sqr_mod(r, p, q);\n" \
"		if (bittest_80(e, b)) r = dup_mod(r, p);\n" \
"	}\n" \
"	return (eq80(r, one) || eq80(r, sub80(p, one)));\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void generate_primes(__global sz_t * restrict const prime_count,\n" \
"	__global uint_64 * restrict const k_vector, __global uint_64_2 * restrict const q_vector,\n" \
"	__global uint_64_2 * restrict const ext_vector, const uint_64 i)\n" \
"{\n" \
"	const sz_t id = (sz_t)get_global_id(0);\n" \
"\n" \
"	const uint_64 j = (i << LN_BLKSZ) | id;\n" \
"	const uint_64 k = 15 * (j / 8) + wheel[j % 8];\n" \
"\n" \
"	uint_80 p; p.lo = (k << G_N) | 1; p.hi = (uint_16)(k >> (64 - G_N));\n" \
"	const uint_80 q = invert(p);\n" \
"	// one is the Montgomery form of 1: 2^80 mod p\n" \
"	const uint_80 one = modinv80(p);\n" \
"\n" \
"	if (prp(p, q, one))\n" \
"	{\n" \
"		const uint prime_index = atomic_inc(prime_count);\n" \
"		k_vector[prime_index] = k;\n" \
"		q_vector[prime_index] = (uint_64_2)(q.lo, q.hi);\n" \
"		ext_vector[prime_index] = (uint_64_2)(one.lo, one.hi);\n" \
"	}\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void init_factors(__global const sz_t * restrict const prime_count,\n" \
"	__global const uint_64 * restrict const k_vector, __global uint_64_2 * restrict const q_vector,\n" \
"	__global uint_64_2 * restrict const ext_vector, __global const int_8 * restrict const kro_vector,\n" \
"	__global uint_64_2 * restrict const c_vector, __global uint_64_2 * restrict const cn_vector)\n" \
"{\n" \
"	const sz_t id = (sz_t)get_global_id(0);\n" \
"	if (id >= *prime_count) return;\n" \
"\n" \
"	const uint_64 k = k_vector[id];\n" \
"	const uint_64_2 ql = q_vector[id], onel = ext_vector[id];\n" \
"	uint_80 q; q.lo = ql.s0; q.hi = (uint_16)(ql.s1);\n" \
"	uint_80 one; one.lo = onel.s0; one.hi = (uint_16)(onel.s1);\n" \
"\n" \
"	uint_80 p; p.lo = (k << G_N) | 1; p.hi = (uint_16)(k >> (64 - G_N));\n" \
"	const uint_80 two = dup_mod(one, p);\n" \
"\n" \
"	// p = 1 (mod 4). If a is odd then (a/p) = (p/a) = ({p mod a}/a)\n" \
"\n" \
"	uint_32 a = 3; uint_80 am = add_mod(two, one, p);\n" \
"	const uint_32 pmod3 = (((uint_32)(k % 3) << G_N) | 1) % 3;\n" \
"	if (pmod3 != 2)\n" \
"	{\n" \
"		a += 2; am = add_mod(am, two, p);\n" \
"		const uint_32 pmod5 = (((uint_32)(k % 5) << G_N) | 1) % 5;\n" \
"		if (kro_vector[256 * ((5 - 3) / 2) + pmod5] >= 0)\n" \
"		{\n" \
"			a += 2; am = add_mod(am, two, p);\n" \
"			const uint_32 pmod7 = (((uint_32)(k % 7) << G_N) | 1) % 7;\n" \
"			if (kro_vector[256 * ((7 - 3) / 2) + pmod7] >= 0)\n" \
"			{\n" \
"				a += 4; am = add_mod(am, dup_mod(two, p), p);\n" \
"				while (a < 256)\n" \
"				{\n" \
"					const uint_32 pmoda = (uint_32)((((k % a) << G_N) | 1) % a);\n" \
"					if (kro_vector[256 * ((a - 3) / 2) + pmoda] < 0) break;\n" \
"					a += 2; am = add_mod(am, two, p);\n" \
"				}\n" \
"				if (a >= 256)\n" \
"				{\n" \
"					c_vector[id] = (uint_64_2)(0, 0);\n" \
"					return;\n" \
"				}\n" \
"			}\n" \
"		}\n" \
"	}\n" \
"\n" \
"	// a^{(p - 1)/2} = -1 <=> a^{k*2^n} = -1. (a^k)^{2*i + 1} are the roots of b^{2^n} + 1 = 0 (mod p)\n" \
"	const uint_80 cm = pow_mod(am, k, p, q);\n" \
"	const uint_80 c = to_int(cm, p, q), c2 = mul_mod(c, cm, p, q);\n" \
"	c_vector[id] = (uint_64_2)(c.lo, c.hi);\n" \
"	ext_vector[id] = (uint_64_2)(c2.lo, c2.hi);\n" \
"	const uint_80 c2p = div80(c2, p);\n" \
"	q_vector[id] = (uint_64_2)(c2p.lo, c2p.hi);\n" \
"	const uint_80 cn = c; //sub80(p, c);\n" \
"	cn_vector[id] = (uint_64_2)(cn.lo, cn.hi);\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void check_factors(__global const sz_t * restrict const prime_count,\n" \
"	__global const uint_64 * restrict const k_vector, __global const uint_64_2 * restrict const q_vector,\n" \
"	__global uint_64_2 * restrict const c_vector, __global const uint_64_2 * restrict const ext_vector,\n" \
"	__global const uint_64_2 * restrict const cn_vector,\n" \
"	__global sz_t * restrict const factor_count, __global uint_64_2 * restrict const factor_vector,\n" \
"	__global sz_t * restrict const error_count, __global uint_64 * restrict const error_vector,\n" \
"	const char last)\n" \
"{\n" \
"	const sz_t id = (sz_t)get_global_id(0);\n" \
"	if (id >= *prime_count) return;\n" \
"\n" \
"	const uint_64_2 cl = c_vector[id];\n" \
"	uint_80 c; c.lo = cl.s0; c.hi = (uint_16)(cl.s1);\n" \
"	if (z80(c)) return;\n" \
"\n" \
"	const uint_64 k = k_vector[id];\n" \
"	uint_80 p; p.lo = (k << G_N) | 1; p.hi = (uint_16)(k >> (64 - G_N));\n" \
"\n" \
"	const uint_64_2 c0sql = ext_vector[id], c0sqpl = q_vector[id];\n" \
"	uint_80 c0sq; c0sq.lo = c0sql.s0; c0sq.hi = (uint_16)(c0sql.s1);\n" \
"	uint_80 c0sqp; c0sqp.lo = c0sqpl.s0; c0sqp.hi = (uint_16)(c0sqpl.s1);\n" \
"\n" \
"	for (size_t l = 0; l < FBLK; ++l)\n" \
"	{\n" \
"		const uint_80 b = (c.lo % 2 == 0) ? c : sub80(p, c);\n" \
"\n" \
"		if ((b.hi == 0) && (b.lo <= 2000000000))\n" \
"		{\n" \
"			const sz_t factor_index = atomic_inc(factor_count);\n" \
"			factor_vector[factor_index] = (uint_64_2)(k, b.lo);\n" \
"		}\n" \
"\n" \
"		c = mul_mod_vs(c, c0sq, c0sqp, p);	// c = a^{(2*i + 1).k}\n" \
"	}\n" \
"\n" \
"	if (last == (char)(0)) c_vector[id] = (uint_64_2)(c.lo, c.hi);\n" \
"	else\n" \
"	{\n" \
"		const uint_64_2 cnegl = cn_vector[id];\n" \
"		uint_80 cn; cn.lo = cnegl.s0; cn.hi = (uint_16)(cnegl.s1);\n" \
"		cn = sub80(p, cn);\n" \
"		if (!eq80(c, cn))\n" \
"		{\n" \
"			const sz_t error_index = atomic_inc(error_count);\n" \
"			error_vector[error_index] = k;\n" \
"		}\n" \
"	}\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void clear(__global sz_t * restrict const count)\n" \
"{\n" \
"	*count = 0;\n" \
"}\n" \
"";
