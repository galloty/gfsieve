/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>

static const char * const src_ocl_sieve64 = \
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
"inline int ilog2(const uint64 x) { return 63 - clz(x); }\n" \
"inline bool bittest(const uint64 x, const int b) { return ((x & ((uint64)(1) << b)) != 0); }\n" \
"\n" \
"inline uint64 add_mod(const uint64 a, const uint64 b, const uint64 p) { return a + b - ((a >= p - b) ? p : 0); }\n" \
"inline uint64 dup_mod(const uint64 a, const uint64 p) { return add_mod(a, a, p); }\n" \
"inline uint64 sub_mod(const uint64 a, const uint64 b, const uint64 p) { return a - b + ((a < b) ? p : 0); }\n" \
"\n" \
"// Montgomery modular multiplication\n" \
"inline uint64 mul_mod(const uint64 a, const uint64 b, const uint64 p, const uint64 q)\n" \
"{\n" \
"	const uint64 ab_l = a * b, ab_h = mul_hi(a, b);\n" \
"	return sub_mod(ab_h, mul_hi(ab_l * q, p), p);\n" \
"}\n" \
"\n" \
"inline uint64 sqr_mod(const uint64 a, const uint64 p, const uint64 q) { return mul_mod(a, a, p, q); }\n" \
"\n" \
"inline uint64 mul_mod_const(const uint64 a, const uint64 b, const uint64 p, const uint64 bq)\n" \
"{\n" \
"	const uint64 ab_h = mul_hi(a, b);\n" \
"	return sub_mod(ab_h, mul_hi(a * bq, p), p);\n" \
"}\n" \
"\n" \
"// Conversion out of Montgomery form\n" \
"inline uint64 Montgomery2int(const uint64 r, const uint64 p, const uint64 q)\n" \
"{\n" \
"	const uint64 mp = mul_hi(r * q, p);\n" \
"	return (mp != 0) ? p - mp : 0;\n" \
"}\n" \
"\n" \
"// p * p_inv = 1 (mod 2^64) (Newton's method)\n" \
"inline uint64 invert(const uint64 p)\n" \
"{\n" \
"	uint64 p_inv = 2 - p;\n" \
"	uint64 prev; do { prev = p_inv; p_inv *= 2 - p * p_inv; } while (p_inv != prev);\n" \
"	return p_inv;\n" \
"}\n" \
"\n" \
"// a^e mod p, left-to-right algorithm\n" \
"inline uint64 pow_mod(const uint64 a, const uint64 e, const uint64 p, const uint64 q)\n" \
"{\n" \
"	const uint64 aq = a * q;\n" \
"	uint64 r = a;\n" \
"	for (int b = ilog2(e) - 1; b >= 0; --b)\n" \
"	{\n" \
"		r = sqr_mod(r, p, q);\n" \
"		if (bittest(e, b)) r = mul_mod_const(r, a, p, aq);\n" \
"	}\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 2^{(p - 1)/2} ?= +/-1 mod p\n" \
"inline bool prp(const uint64 p, const uint64 q, const uint64 one)\n" \
"{\n" \
"	const uint64 e = (p - 1) / 2;\n" \
"	int b = ilog2(e) - 1;\n" \
"	uint64 r = dup_mod(one, p);	// 2 = 1 + 1\n" \
"	r = dup_mod(r, p);			// 2 * 2 = 2 + 2\n" \
"	if (bittest(e, b)) r = dup_mod(r, p);\n" \
"	for (--b; b >= 0; --b)\n" \
"	{\n" \
"		r = sqr_mod(r, p, q);\n" \
"		if (bittest(e, b)) r = dup_mod(r, p);\n" \
"	}\n" \
"	return ((r == one) || (r == p - one));\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void generate_primes(__global sz_t * restrict const prime_count,\n" \
"	__global uint64 * restrict const k_vector, __global uint64 * restrict const q_vector,\n" \
"	__global uint64 * restrict const ext_vector, const uint64 i)\n" \
"{\n" \
"	const sz_t id = (sz_t)get_global_id(0);\n" \
"\n" \
"	const uint64 j = (i << LN_BLKSZ) | id;\n" \
"	const uint64 k = 15 * (j / 8) + wheel[j % 8];\n" \
"\n" \
"	// one is the Montgomery form of 1: 2^64 mod p = (2^64 - p) mod p\n" \
"	const uint64 p = (k << G_N) | 1, q = invert(p), one = (-p) % p;\n" \
"	if (prp(p, q, one))\n" \
"	{\n" \
"		const uint prime_index = atomic_inc(prime_count);\n" \
"		k_vector[prime_index] = k;\n" \
"		q_vector[prime_index] = q;\n" \
"		ext_vector[prime_index] = one;\n" \
"	}\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void init_factors(__global const sz_t * restrict const prime_count,\n" \
"	__global const uint64 * restrict const k_vector, __global const uint64 * restrict const q_vector,\n" \
"	__global uint64 * restrict const ext_vector, __global const int_8 * restrict const kro_vector,\n" \
"	__global uint64 * restrict const c_vector)\n" \
"{\n" \
"	const sz_t id = (sz_t)get_global_id(0);\n" \
"	if (id >= *prime_count) return;\n" \
"\n" \
"	const uint64 k = k_vector[id], q = q_vector[id], one = ext_vector[id];\n" \
"\n" \
"	const uint64 p = (k << G_N) | 1;\n" \
"	const uint64 two = dup_mod(one, p);\n" \
"\n" \
"	// p = 1 (mod 4). If a is odd then (a/p) = (p/a) = ({p mod a}/a)\n" \
"\n" \
"	uint32 a = 3; uint64 am = add_mod(two, one, p);\n" \
"	if (p % 3 != 2)\n" \
"	{\n" \
"		a += 2; am = add_mod(am, two, p);\n" \
"		if (kro_vector[256 * ((5 - 3) / 2) + (uint32)(p % 5)] >= 0)\n" \
"		{\n" \
"			a += 2; am = add_mod(am, two, p);\n" \
"			if (kro_vector[256 * ((7 - 3) / 2) + (uint32)(p % 7)] >= 0)\n" \
"			{\n" \
"				a += 4; am = add_mod(am, dup_mod(two, p), p);\n" \
"				while (a < 256)\n" \
"				{\n" \
"					if (kro_vector[256 * ((a - 3) / 2) + (uint32)(p % a)] < 0) break;\n" \
"					a += 2; am = add_mod(am, two, p);\n" \
"				}\n" \
"				if (a >= 256)\n" \
"				{\n" \
"					c_vector[id] = 0;\n" \
"					return;\n" \
"				}\n" \
"			}\n" \
"		}\n" \
"	}\n" \
"\n" \
"	// a^{(p - 1)/2} = -1 <=> a^{k*2^n} = -1. (a^k)^{2*i + 1} are the roots of b^{2^n} + 1 = 0 (mod p)\n" \
"	const uint64 cm = pow_mod(am, k, p, q);\n" \
"	c_vector[id] = Montgomery2int(cm, p, q);\n" \
"	ext_vector[id] = sqr_mod(cm, p, q);\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void check_factors(__global const sz_t * restrict const prime_count,\n" \
"	__global const uint64 * restrict const k_vector, __global const uint64 * restrict const q_vector,\n" \
"	__global uint64 * restrict const c_vector, __global const uint64 * restrict const ext_vector,\n" \
"	__global sz_t * restrict const factor_count, __global uint64_2 * restrict const factor_vector)\n" \
"{\n" \
"	const sz_t id = (sz_t)get_global_id(0);\n" \
"	if (id >= *prime_count) return;\n" \
"\n" \
"	uint64 c = c_vector[id];\n" \
"	if (c == 0) return;\n" \
"\n" \
"	const uint64 k = k_vector[id], p = (k << G_N) | 1, q = q_vector[id];\n" \
"	const uint64 c0sq = ext_vector[id], c0sqxq = c0sq * q;\n" \
"\n" \
"	for (size_t l = 0; l < FBLK; ++l)\n" \
"	{\n" \
"		const uint64 b = (c % 2 == 0) ? c : p - c;\n" \
"\n" \
"		if (b <= 2000000000)\n" \
"		{\n" \
"			const sz_t factor_index = atomic_inc(factor_count);\n" \
"			factor_vector[factor_index] = (uint64_2)(k, b);\n" \
"\n" \
"		}\n" \
"\n" \
"		c = mul_mod_const(c, c0sq, p, c0sqxq);	// c = a^{(2*i + 1).k}\n" \
"	}\n" \
"\n" \
"	c_vector[id] = c;\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void clear(__global sz_t * restrict const count)\n" \
"{\n" \
"	*count = 0;\n" \
"}\n" \
"";
