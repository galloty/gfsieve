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
"typedef uint	uint32;\n" \
"typedef ulong	uint64;\n" \
"\n" \
"inline int ilog2(const uint64 x) { return 63 - clz(x); }\n" \
"\n" \
"inline bool bittest(const uint64 x, const int b) { return ((x & ((uint64)(1) << b)) != 0); }\n" \
"\n" \
"inline uint64 add_mod(const uint64 a, const uint64 b, const uint64 p)\n" \
"{\n" \
"	const uint64 c = (a >= p - b) ? p : 0;\n" \
"	return a + b - c;\n" \
"}\n" \
"\n" \
"inline uint64 sub_mod(const uint64 a, const uint64 b, const uint64 p)\n" \
"{\n" \
"	const uint64 c = (a < b) ? p : 0;\n" \
"	return a - b + c;\n" \
"}\n" \
"\n" \
"inline uint64 mul_mod(const uint64 a, const uint64 b, const uint64 p, const uint64 q)\n" \
"{\n" \
"	const uint64 ab_l = a * b, ab_h = mul_hi(a, b);\n" \
"	const uint64 mp = mul_hi(ab_l * q, p), r = ab_h - mp;\n" \
"	return (ab_h < mp) ? r + p : r;\n" \
"}\n" \
"\n" \
"inline uint64 mul_mod_const(const uint64 a, const uint64 b, const uint64 p, const uint64 q, const uint64 bq)\n" \
"{\n" \
"	const uint64 ab_h = mul_hi(a, b);\n" \
"	const uint64 mp = mul_hi(a * bq, p), r = ab_h - mp;\n" \
"	return (ab_h < mp) ? r + p : r;\n" \
"}\n" \
"\n" \
"inline uint64 toInt(const uint64 r, const uint64 p, const uint64 q)\n" \
"{\n" \
"	const uint64 mp = mul_hi(r * q, p);\n" \
"	return (mp != 0) ? p - mp : 0;\n" \
"}\n" \
"\n" \
"// p * p_inv = 1 (mod 2^64) (Newton's method)\n" \
"inline uint64 invert(const uint64 p)\n" \
"{\n" \
"	uint64 p_inv = 1, prev = 0;\n" \
"	while (p_inv != prev) { prev = p_inv; p_inv *= 2 - p * p_inv; }\n" \
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
"		r = mul_mod(r, r, p, q);\n" \
"		if (bittest(e, b)) r = mul_mod_const(r, a, p, q, aq);\n" \
"	}\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 2^{(p - 1)/2} ?= +/-1 mod p\n" \
"inline bool prp(const uint64 p, const uint64 q, const uint64 one)\n" \
"{\n" \
"	const uint64 e = (p - 1) / 2;\n" \
"	int b = ilog2(e) - 1;\n" \
"	uint64 r = add_mod(one, one, p); r = add_mod(r, r, p);\n" \
"	if (bittest(e, b)) r = add_mod(r, r, p);\n" \
"	for (--b; b >= 0; --b)\n" \
"	{\n" \
"		r = mul_mod(r, r, p, q);\n" \
"		if (bittest(e, b)) r = add_mod(r, r, p);\n" \
"	}\n" \
"	return ((r == one) || (r == p - one));\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void generate_primes(__global uint * restrict const prime_count, __global ulong * restrict const k_vector,\n" \
"	__global ulong2 * restrict const k_ext_vector, const ulong i)\n" \
"{\n" \
"	const uint64 k = (i << log2GlobalWorkSize) | get_global_id(0);\n" \
"\n" \
"	const uint64 p = (k << g_n) | 1, q = invert(p), one = (-p) % p;\n" \
"	if (prp(p, q, one))\n" \
"	{\n" \
"		const uint prime_index = atomic_inc(prime_count);\n" \
"		k_vector[prime_index] = k;\n" \
"		k_ext_vector[prime_index] = (ulong2)(q, one);\n" \
"	}\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void init_factors(__global const uint * restrict const prime_count, __global const ulong * restrict const k_vector,\n" \
"	__global const ulong2 * restrict const k_ext_vector, __global const char * restrict const _kro_vector,\n" \
"	__global ulong2 * restrict const c_vector, __global ulong2 * restrict const a2k_vector)\n" \
"{\n" \
"	const size_t i = get_global_id(0);\n" \
"	if (i >= *prime_count) return;\n" \
"\n" \
"	const uint64 k = k_vector[i];\n" \
"	const ulong2 k_ext = k_ext_vector[i];\n" \
"\n" \
"	const uint64 p = (k << g_n) | 1, q = k_ext.s0, one = k_ext.s1;\n" \
"	const uint64 two = add_mod(one, one, p);\n" \
"\n" \
"	// p = 1 (mod 4). If a is odd then (a/p) = (p/a) = ({p mod a}/a)\n" \
"\n" \
"	uint32 a = 3; uint64 am = add_mod(two, one, p);\n" \
"	if (p % 3 != 2)\n" \
"	{\n" \
"		a += 2; am = add_mod(am, two, p);\n" \
"		if (_kro_vector[256 * ((5 - 3) / 2) + (uint32)(p % 5)] >= 0)\n" \
"		{\n" \
"			a += 2; am = add_mod(am, two, p);\n" \
"			if (_kro_vector[256 * ((7 - 3) / 2) + (uint32)(p % 7)] >= 0)\n" \
"			{\n" \
"				a += 4; am = add_mod(am, two, p); am = add_mod(am, two, p);\n" \
"				while (a < 256)\n" \
"				{\n" \
"					if (_kro_vector[256 * ((a - 3) / 2) + (uint32)(p % a)] < 0) break;\n" \
"					a += 2; am = add_mod(am, two, p);\n" \
"				}\n" \
"				if (a >= 256)\n" \
"				{\n" \
"					c_vector[i] = (ulong2)(0, 0);\n" \
"					return;\n" \
"				}\n" \
"			}\n" \
"		}\n" \
"	}\n" \
"\n" \
"	const uint64 cm = pow_mod(am, k, p, q);\n" \
"\n" \
"	const uint64 c = toInt(cm, p, q);\n" \
"	c_vector[i] = (ulong2)(c, 0);\n" \
"\n" \
"	const uint64 a2km = mul_mod(cm, cm, p, q);\n" \
"	a2k_vector[i] = (ulong2)(a2km, 0);\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void check_factors(__global const uint * restrict const prime_count,\n" \
"	__global const ulong * restrict const k_vector, __global const ulong2 * restrict const k_ext_vector,\n" \
"	__global ulong2 * restrict const c_vector, __global const ulong2 * restrict const a2k_vector,\n" \
"	__global uint * restrict const factor_count, __global ulong2 * restrict const factor)\n" \
"{\n" \
"	const size_t i = get_global_id(0);\n" \
"	if (i >= *prime_count) return;\n" \
"\n" \
"	const uint64 p = (k_vector[i] << g_n) | 1, q = k_ext_vector[i].s0;\n" \
"\n" \
"	uint64 c = c_vector[i].s0;\n" \
"	if (c == 0) return;\n" \
"	const uint64 a2k = a2k_vector[i].s0, a2kq = a2k * q;\n" \
"\n" \
"	for (size_t i = 0; i < factors_loop; ++i)\n" \
"	{\n" \
"		if (c % 2 != 0) c = p - c;\n" \
"		if (c <= 2000000000)\n" \
"		{\n" \
"			const uint factor_index = atomic_inc(factor_count);\n" \
"			factor[factor_index] = (ulong2)(p, c << 32);\n" \
"		}\n" \
"		c = mul_mod_const(c, a2k, p, q, a2kq);		// c = a^{(2*i + 1).k}\n" \
"	}\n" \
"\n" \
"	c_vector[i].s0 = c;\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void clear_primes(__global uint * restrict const prime_count)\n" \
"{\n" \
"	*prime_count = 0;\n" \
"}\n" \
"";
