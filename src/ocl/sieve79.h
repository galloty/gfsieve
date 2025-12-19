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
"#define sz_t		uint\n" \
"#ifndef uint_8\n" \
"#define uint_8		uchar\n" \
"#endif\n" \
"#define int_8		char\n" \
"#define uint16		ushort\n" \
"#define uint32		uint\n" \
"#define uint64		ulong\n" \
"#define uint64_2	ulong2\n" \
"\n" \
"#if !defined(G_N)\n" \
"#define G_N			16\n" \
"#define LN_BLKSZ	22\n" \
"#define FBLK		1024\n" \
"\n" \
"__constant uint_8 wheel[8] = { 0, 1, 3, 6, 7, 10, 12, 13 };\n" \
"#endif\n" \
"\n" \
"typedef struct _uint80\n" \
"{\n" \
"	uint64 s0;\n" \
"	uint16 s1;\n" \
"} uint80;\n" \
"\n" \
"inline int ilog2_64(const uint64 x) { return 63 - clz(x); }\n" \
"inline int ilog2_80(const uint80 x) { return (x.s1 != 0) ? (64 + 31 - clz((uint32)x.s1)) : ilog2_64(x.s0); }\n" \
"inline bool bittest_64(const uint64 x, const int b) { return ((x & ((uint64)(1) << b)) != 0); }\n" \
"inline bool bittest_80(const uint80 x, const int b) { return (b >= 64) ? ((x.s1 & ((uint16)(1) << (b - 64))) != 0) : bittest_64(x.s0, b); }\n" \
"\n" \
"inline uint32 lo32(const uint64 x) { return (uint32)(x); }\n" \
"inline uint32 hi32(const uint64 x) { return (uint32)(x >> 32); }\n" \
"\n" \
"inline bool eq80(const uint80 x, const uint80 y)	// x is equal to y\n" \
"{\n" \
"	return ((x.s1 == y.s1) && (x.s0 == y.s0));\n" \
"}\n" \
"\n" \
"inline bool ge80(const uint80 x, const uint80 y)	// x is greater than or equal to y\n" \
"{\n" \
"	if (x.s1 > y.s1) return true;\n" \
"	if (x.s1 < y.s1) return false;\n" \
"	return (x.s0 >= y.s0);\n" \
"}\n" \
"\n" \
"inline bool z80(const uint80 x)	// x is equal to 0\n" \
"{\n" \
"	return ((x.s1 == 0) && (x.s0 == 0));\n" \
"}\n" \
"\n" \
"inline uint80 add80(const uint80 x, const uint80 y)\n" \
"{\n" \
"	uint80 r;\n" \
"#ifdef PTX_ASM\n" \
"	const uint32 xs1 = x.s1, ys1 = y.s1;\n" \
"	uint32 rs1;\n" \
"	asm volatile (\"add.cc.u64 %0, %1, %2;\" : \"=l\" (r.s0) : \"l\" (x.s0), \"l\" (y.s0));\n" \
"	asm volatile (\"addc.u32 %0, %1, %2;\" : \"=r\" (rs1) : \"r\" (xs1), \"r\" (ys1));\n" \
"	r.s1 = (uint16)(rs1);\n" \
"#else\n" \
"	r.s0 = x.s0 + y.s0; r.s1 = x.s1 + y.s1 + ((r.s0 < x.s0) ? 1 : 0);\n" \
"#endif\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint80 sub80(const uint80 x, const uint80 y)\n" \
"{\n" \
"	uint80 r;\n" \
"#ifdef PTX_ASM\n" \
"	const uint32 xs1 = x.s1, ys1 = y.s1;\n" \
"	uint32 rs1;\n" \
"	asm volatile (\"sub.cc.u64 %0, %1, %2;\" : \"=l\" (r.s0) : \"l\" (x.s0), \"l\" (y.s0));\n" \
"	asm volatile (\"subc.u32 %0, %1, %2;\" : \"=r\" (rs1) : \"r\" (xs1), \"r\" (ys1));\n" \
"	r.s1 = (uint16)(rs1);\n" \
"#else\n" \
"	 r.s0 = x.s0 - y.s0; r.s1 = x.s1 - y.s1 - ((x.s0 < y.s0) ? 1 : 0);\n" \
"#endif\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint80 neg80(const uint80 x)\n" \
"{\n" \
"	uint80 r;\n" \
"#ifdef PTX_ASM\n" \
"	const uint32 xs1 = x.s1;\n" \
"	uint32 rs1;\n" \
"	asm volatile (\"sub.cc.u64 %0, 0, %1;\" : \"=l\" (r.s0) : \"l\" (x.s0));\n" \
"	asm volatile (\"subc.u32 %0, 0, %1;\" : \"=r\" (rs1) : \"r\" (xs1));\n" \
"	r.s1 = (uint16)(rs1);\n" \
"#else\n" \
"	r.s0 = -x.s0; r.s1 = -x.s1 - ((x.s0 != 0) ? 1 : 0);\n" \
"#endif\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 0 < s < 64\n" \
"inline uint80 shl80(const uint80 x, const int s)\n" \
"{\n" \
"	const uint16 rs1 = (s < 16) ? (x.s1 << s) : 0;\n" \
"	uint80 r; r.s0 = x.s0 << s; r.s1 = rs1 | (uint16)(x.s0 >> (64 - s));\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 0 < s < 64\n" \
"inline uint80 shr80(const uint80 x, const int s)\n" \
"{\n" \
"	const uint16 rs1 = (s < 16) ? (x.s1 >> s) : 0;\n" \
"	uint80 r; r.s0 = (x.s0 >> s) | ((uint64)(x.s1) << (64 - s)); r.s1 = rs1;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"typedef struct _uint96\n" \
"{\n" \
"	uint64 s0;\n" \
"	uint32 s1;\n" \
"} uint96;\n" \
"\n" \
"inline uint96 madd96(const uint96 z, const uint64 x, const uint32 y)\n" \
"{\n" \
"	uint96 r;\n" \
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
"	r.s0 = z.s0 + x * y; r.s1 = z.s1 + (uint32)(mul_hi(x, (uint64)(y))) + ((r.s0 < z.s0) ? 1 : 0);\n" \
"#endif\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline void mul80_wide(const uint80 x, const uint80 y, uint80 * const lo, uint80 * const hi)\n" \
"{\n" \
"	lo->s0 = x.s0 * y.s0;\n" \
"	uint96 t; t.s0 = mul_hi(x.s0, y.s0); t.s1 = x.s1 * (uint32)(y.s1);\n" \
"	const uint96 r = madd96(madd96(t, x.s0, y.s1), y.s0, x.s1);\n" \
"	lo->s1 = (uint16)(r.s0); hi->s0 = (r.s0 >> 16) | ((uint64)(r.s1) << 48); hi->s1 = (uint16)(r.s1 >> 16);\n" \
"}\n" \
"\n" \
"inline void sqr80_wide(const uint80 x, uint80 * const lo, uint80 * const hi)\n" \
"{\n" \
"	lo->s0 = x.s0 * x.s0;\n" \
"	uint96 t; t.s0 = mul_hi(x.s0, x.s0); t.s1 = x.s1 * (uint32)(x.s1);\n" \
"	const uint96 r = madd96(t, x.s0, (uint32)(x.s1) << 1);\n" \
"	lo->s1 = (uint16)(r.s0); hi->s0 = (r.s0 >> 16) | ((uint64)(r.s1) << 48); hi->s1 = (uint16)(r.s1 >> 16);\n" \
"}\n" \
"\n" \
"inline uint80 mul80(const uint80 x, const uint80 y)\n" \
"{\n" \
"	const uint32 a0 = (uint32)(x.s0), a1 = (uint32)(x.s0 >> 32), a2 = x.s1;\n" \
"	const uint32 b0 = (uint32)(y.s0), b1 = (uint32)(y.s0 >> 32), b2 = y.s1;\n" \
"	uint32 c0 = a0 * b0, c1 = mul_hi(a0, b0), c2 = a1 * b1;\n" \
"#ifdef PTX_ASM\n" \
"	asm volatile (\"mad.lo.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a0), \"r\" (b2), \"r\" (c2));\n" \
"	asm volatile (\"mad.lo.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a2), \"r\" (b0), \"r\" (c2));\n" \
"	asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c1) : \"r\" (a0), \"r\" (b1), \"r\" (c1));\n" \
"	asm volatile (\"madc.hi.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a0), \"r\" (b1), \"r\" (c2));\n" \
"	asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c1) : \"r\" (a1), \"r\" (b0), \"r\" (c1));\n" \
"	asm volatile (\"madc.hi.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a1), \"r\" (b0), \"r\" (c2));\n" \
"#else\n" \
"	c2 += a0 * b2 + a2 * b0;\n" \
"	const uint64 c12 = upsample(c2, c1) + a0 * (uint64)(b1) + a1 * (uint64)(b0);\n" \
"	c1 = lo32(c12); c2 = hi32(c12);\n" \
"#endif\n" \
"	uint80 r; r.s0 = upsample(c1, c0); r.s1 = (uint16)(c2);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint80 mul80_hi(const uint80 x, const uint80 y)\n" \
"{\n" \
"	uint96 t; t.s0 = mul_hi(x.s0, y.s0); t.s1 = x.s1 * (uint32)(y.s1);\n" \
"	const uint96 r96 = madd96(madd96(t, x.s0, y.s1), y.s0, x.s1);\n" \
"	uint80 r; r.s0 = (r96.s0 >> 16) | ((uint64)(r96.s1) << 48); r.s1 = (uint16)(r96.s1 >> 16);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"typedef struct _uint160\n" \
"{\n" \
"	uint64 s0, s1;\n" \
"	uint32 s2;\n" \
"} uint160;\n" \
"\n" \
"inline bool ge160(const uint160 x, const uint160 y)	// x is greater than or equal to y\n" \
"{\n" \
"	if (x.s2 > y.s2) return true;\n" \
"	if (x.s2 < y.s2) return false;\n" \
"	if (x.s1 > y.s1) return true;\n" \
"	if (x.s1 < y.s1) return false;\n" \
"	return (x.s0 >= y.s0);\n" \
"}\n" \
"\n" \
"inline uint160 sub160(const uint160 x, const uint160 y)\n" \
"{\n" \
"	uint160 r;\n" \
"#ifdef PTX_ASM\n" \
"	asm volatile (\"sub.cc.u64 %0, %1, %2;\" : \"=l\" (r.s0) : \"l\" (x.s0), \"l\" (y.s0));\n" \
"	asm volatile (\"subc.cc.u64 %0, %1, %2;\" : \"=l\" (r.s1) : \"l\" (x.s1), \"l\" (y.s1));\n" \
"	asm volatile (\"subc.u32 %0, %1, %2;\" : \"=r\" (r.s2) : \"r\" (x.s2), \"r\" (y.s2));\n" \
"#else\n" \
"	const uint32 c0 = (x.s0 < y.s0) ? 1 : 0, c1 = (x.s1 < y.s1) ? 1 : 0;\n" \
"	r.s0 = x.s0 - y.s0; r.s1 = x.s1 - y.s1; r.s2 = x.s2 - y.s2;\n" \
"	const uint32 c2 = (r.s1 < c0) ? 1 : 0;\n" \
"	r.s1 -= c0; r.s2 -= c1 + c2;\n" \
"#endif\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 0 < s < 64\n" \
"inline uint160 shl160(const uint160 x, const int s)\n" \
"{\n" \
"	const uint32 rs2 = (s < 32) ? (x.s2 << s) : 0;\n" \
"	uint160 r; r.s0 = x.s0 << s; r.s1 = (x.s1 << s) | (x.s0 >> (64 - s)); r.s2 = rs2 | (uint32)(x.s1 >> (64 - s));\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 0 < s < 64\n" \
"inline uint160 shr160(const uint160 x, const int s)\n" \
"{\n" \
"	const uint32 rs2 = (s < 32) ? (x.s2 >> s) : 0;\n" \
"	uint160 r; r.s0 = (x.s0 >> s) | (x.s1 << (64 - s)); r.s1 = (x.s1 >> s) | ((uint64)(x.s2) << (64 - s)); r.s2 = rs2;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 2^80 * x / y, x < y, y is not a power of 2\n" \
"inline uint80 div80(const uint80 x, const uint80 y)\n" \
"{\n" \
"	if (z80(x)) return x;\n" \
"\n" \
"	// 2^b * y <= 2^80 * x < 2^{b+1} * y, 0 <= b < 80\n" \
"	int b = 80 + ilog2_80(x) - ilog2_80(y) - 1;\n" \
"\n" \
"	// m = 2^b * y\n" \
"	const int b_64 = (b < 64) ? 0 : 1, b64 = b % 64;\n" \
"	uint160 m; m.s0 = (b_64 == 0) ? y.s0 : 0; m.s1 = (b_64 == 0) ? (uint64)(y.s1) : y.s0; m.s2 = (b_64 == 0) ? 0 : (uint32)(y.s1);\n" \
"	if (b64 > 0) m = shl160(m, b64);\n" \
"\n" \
"	// z = 2^80 * x\n" \
"	uint160 z; z.s0 = 0; z.s1 = x.s0 << 16; z.s2 = ((uint32)(x.s1) << 16) | (uint32)(x.s0 >> 48);\n" \
"	// z < 2*m ?\n" \
"	uint160 t = shl160(m, 1);\n" \
"	if (ge160(z, t)) { ++b; m = t; }\n" \
"\n" \
"	uint80 r; r.s0 = 0; r.s1 = 0;\n" \
"	for (int j = 0; j <= b; ++j)\n" \
"	{\n" \
"		r = shl80(r, 1);\n" \
"		if (ge160(z, m)) { z = sub160(z, m); r.s0 |= 1; }\n" \
"		m = shr160(m, 1);\n" \
"	}\n" \
"\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 2^80 mod x, x is not a power of 2, 2^63 < x < 2^79\n" \
"inline uint80 modinv80(const uint80 x)\n" \
"{\n" \
"	// 2^b * x < 2^80 < 2^{b+1} * x, 1 <= b <= 16\n" \
"	const int b = 80 - ilog2_80(x) - 1;\n" \
"\n" \
"	// m = 2^b * x\n" \
"	uint80 m = shl80(x, b);\n" \
"	// r = 2^80 - m\n" \
"	uint80 r = neg80(m);\n" \
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
"inline uint80 add_mod(const uint80 a, const uint80 b, const uint80 p)\n" \
"{\n" \
"	uint80 r = add80(a, b);\n" \
"	if (ge80(r, p)) r = sub80(r, p);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint80 dup_mod(const uint80 a, const uint80 p)\n" \
"{\n" \
"	uint80 r = shl80(a, 1);\n" \
"	if (ge80(r, p)) r = sub80(r, p);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint80 sub_mod(const uint80 a, const uint80 b, const uint80 p)\n" \
"{\n" \
"	uint80 r = sub80(a, b);\n" \
"	if (!ge80(a, b)) r = add80(r, p);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// Montgomery modular multiplication\n" \
"inline uint80 mul_mod(const uint80 a, const uint80 b, const uint80 p, const uint80 q)\n" \
"{\n" \
"	uint80 ab_l, ab_h; mul80_wide(a, b, &ab_l, &ab_h);\n" \
"	return sub_mod(ab_h, mul80_hi(mul80(ab_l, q), p), p);\n" \
"}\n" \
"\n" \
"inline uint80 sqr_mod(const uint80 a, const uint80 p, const uint80 q)\n" \
"{\n" \
"	uint80 a2_l, a2_h; sqr80_wide(a, &a2_l, &a2_h);\n" \
"	return sub_mod(a2_h, mul80_hi(mul80(a2_l, q), p), p);\n" \
"}\n" \
"\n" \
"// Victor Shoup’s modular multiplication\n" \
"inline uint80 mul_mod_vs(const uint80 a, const uint80 b, const uint80 bp, const uint80 p)\n" \
"{\n" \
"	const uint80 abp_h = mul80_hi(a, bp);\n" \
"	const uint80 r = sub80(mul80(a, b), mul80(abp_h, p));\n" \
"	return ge80(r, p) ? sub80(r, p) : r;\n" \
"}\n" \
"\n" \
"inline uint80 to_int(const uint80 a, const uint80 p, const uint80 q)\n" \
"{\n" \
"	const uint80 mp = mul80_hi(mul80(a, q), p);\n" \
"	return !z80(mp) ? sub80(p, mp) : mp;\n" \
"}\n" \
"\n" \
"// p * p_inv = 1 (mod 2^80) (Newton's method)\n" \
"inline uint80 invert(const uint80 p)\n" \
"{\n" \
"	uint80 two; two.s0 = 2; two.s1 = 0;\n" \
"	uint80 p_inv = sub80(two, p);\n" \
"	uint80 prev; do { prev = p_inv; p_inv = mul80(p_inv, sub80(two, mul80(p, p_inv))); } while (!eq80(p_inv, prev));\n" \
"	return p_inv;\n" \
"}\n" \
"\n" \
"// a^e mod p, left-to-right algorithm\n" \
"inline uint80 pow_mod(const uint80 a, const uint64 e, const uint80 p, const uint80 q)\n" \
"{\n" \
"	uint80 r = a;\n" \
"	for (int b = ilog2_64(e) - 1; b >= 0; --b)\n" \
"	{\n" \
"		r = sqr_mod(r, p, q);\n" \
"		if (bittest_64(e, b)) r = mul_mod(r, a, p, q);\n" \
"	}\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 2^{(p - 1)/2} ?= +/-1 mod p\n" \
"inline bool prp(const uint80 p, const uint80 q, const uint80 one)\n" \
"{\n" \
"	uint80 e; e.s0 = (p.s0 >> 1) | ((uint64)(p.s1) << 63); e.s1 = p.s1 >> 1;\n" \
"	int b = ilog2_80(e) - 1;\n" \
"	uint80 r = dup_mod(one, p);	// 2 = 1 + 1\n" \
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
"	__global uint64 * restrict const k_vector, __global uint64_2 * restrict const q_vector,\n" \
"	__global uint64_2 * restrict const ext_vector, const uint64 i)\n" \
"{\n" \
"	const sz_t id = (sz_t)get_global_id(0);\n" \
"\n" \
"	const uint64 j = (i << LN_BLKSZ) | id;\n" \
"	const uint64 k = 15 * (j / 8) + wheel[j % 8];\n" \
"\n" \
"	uint80 p; p.s0 = (k << G_N) | 1; p.s1 = (uint16)(k >> (64 - G_N));\n" \
"	const uint80 q = invert(p);\n" \
"	// one is the Montgomery form of 1: 2^80 mod p\n" \
"	const uint80 one = modinv80(p);\n" \
"\n" \
"	if (prp(p, q, one))\n" \
"	{\n" \
"		const uint prime_index = atomic_inc(prime_count);\n" \
"		k_vector[prime_index] = k;\n" \
"		q_vector[prime_index] = (uint64_2)(q.s0, q.s1);\n" \
"		ext_vector[prime_index] = (uint64_2)(one.s0, one.s1);\n" \
"	}\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void init_factors(__global const sz_t * restrict const prime_count,\n" \
"	__global const uint64 * restrict const k_vector, __global uint64_2 * restrict const q_vector,\n" \
"	__global uint64_2 * restrict const ext_vector, __global const int_8 * restrict const kro_vector,\n" \
"	__global uint64_2 * restrict const c_vector)\n" \
"{\n" \
"	const sz_t id = (sz_t)get_global_id(0);\n" \
"	if (id >= *prime_count) return;\n" \
"\n" \
"	const uint64 k = k_vector[id];\n" \
"	const uint64_2 ql = q_vector[id], onel = ext_vector[id];\n" \
"	uint80 q; q.s0 = ql.s0; q.s1 = (uint16)(ql.s1);\n" \
"	uint80 one; one.s0 = onel.s0; one.s1 = (uint16)(onel.s1);\n" \
"\n" \
"	uint80 p; p.s0 = (k << G_N) | 1; p.s1 = (uint16)(k >> (64 - G_N));\n" \
"	const uint80 two = dup_mod(one, p);\n" \
"\n" \
"	// p = 1 (mod 4). If a is odd then (a/p) = (p/a) = ({p mod a}/a)\n" \
"\n" \
"	uint32 a = 3; uint80 am = add_mod(two, one, p);\n" \
"	const uint32 pmod3 = (((uint32)(k % 3) << G_N) | 1) % 3;\n" \
"	if (pmod3 != 2)\n" \
"	{\n" \
"		a += 2; am = add_mod(am, two, p);\n" \
"		const uint32 pmod5 = (((uint32)(k % 5) << G_N) | 1) % 5;\n" \
"		if (kro_vector[256 * ((5 - 3) / 2) + pmod5] >= 0)\n" \
"		{\n" \
"			a += 2; am = add_mod(am, two, p);\n" \
"			const uint32 pmod7 = (((uint32)(k % 7) << G_N) | 1) % 7;\n" \
"			if (kro_vector[256 * ((7 - 3) / 2) + pmod7] >= 0)\n" \
"			{\n" \
"				a += 4; am = add_mod(am, dup_mod(two, p), p);\n" \
"				while (a < 256)\n" \
"				{\n" \
"					const uint32 pmoda = (uint32)((((k % a) << G_N) | 1) % a);\n" \
"					if (kro_vector[256 * ((a - 3) / 2) + pmoda] < 0) break;\n" \
"					a += 2; am = add_mod(am, two, p);\n" \
"				}\n" \
"				if (a >= 256)\n" \
"				{\n" \
"					c_vector[id] = (uint64_2)(0, 0);\n" \
"					return;\n" \
"				}\n" \
"			}\n" \
"		}\n" \
"	}\n" \
"\n" \
"	// a^{(p - 1)/2} = -1 <=> a^{k*2^n} = -1. (a^k)^{2*i + 1} are the roots of b^{2^n} + 1 = 0 (mod p)\n" \
"	const uint80 cm = pow_mod(am, k, p, q);\n" \
"	const uint80 c = to_int(cm, p, q), c2 = mul_mod(c, cm, p, q);\n" \
"	c_vector[id] = (uint64_2)(c.s0, c.s1);\n" \
"	ext_vector[id] = (uint64_2)(c2.s0, c2.s1);\n" \
"	const uint80 c2p = div80(c2, p);\n" \
"	q_vector[id] = (uint64_2)(c2p.s0, c2p.s1);\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void check_factors(__global const sz_t * restrict const prime_count,\n" \
"	__global const uint64 * restrict const k_vector, __global const uint64_2 * restrict const q_vector,\n" \
"	__global uint64_2 * restrict const c_vector, __global const uint64_2 * restrict const ext_vector,\n" \
"	__global sz_t * restrict const factor_count, __global uint64_2 * restrict const factor_vector)\n" \
"{\n" \
"	const sz_t id = (sz_t)get_global_id(0);\n" \
"	if (id >= *prime_count) return;\n" \
"\n" \
"	const uint64_2 cl = c_vector[id];\n" \
"	uint80 c; c.s0 = cl.s0; c.s1 = (uint16)(cl.s1);\n" \
"	if (z80(c)) return;\n" \
"\n" \
"	const uint64 k = k_vector[id];\n" \
"	uint80 p; p.s0 = (k << G_N) | 1; p.s1 = (uint16)(k >> (64 - G_N));\n" \
"\n" \
"	const uint64_2 c0sql = ext_vector[id], c0sqpl = q_vector[id];\n" \
"	uint80 c0sq; c0sq.s0 = c0sql.s0; c0sq.s1 = (uint16)(c0sql.s1);\n" \
"	uint80 c0sqp; c0sqp.s0 = c0sqpl.s0; c0sqp.s1 = (uint16)(c0sqpl.s1);\n" \
"\n" \
"	for (size_t l = 0; l < FBLK; ++l)\n" \
"	{\n" \
"		const uint80 b = (c.s0 % 2 == 0) ? c : sub80(p, c);\n" \
"\n" \
"		if ((b.s1 == 0) && (b.s0 <= 2000000000))\n" \
"		{\n" \
"			const sz_t factor_index = atomic_inc(factor_count);\n" \
"			factor_vector[factor_index] = (uint64_2)(k, b.s0);\n" \
"		}\n" \
"\n" \
"		c = mul_mod_vs(c, c0sq, c0sqp, p);	// c = a^{(2*i + 1).k}\n" \
"	}\n" \
"\n" \
"	c_vector[id] = (uint64_2)(c.s0, c.s1);\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void clear(__global sz_t * restrict const count)\n" \
"{\n" \
"	*count = 0;\n" \
"}\n" \
"";
