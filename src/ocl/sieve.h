/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>

static const char * const src_ocl_sieve = \
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
"typedef struct _uint96\n" \
"{\n" \
"	uint64 s0;\n" \
"	uint32 s1;\n" \
"} uint96;\n" \
"\n" \
"typedef struct _uint128\n" \
"{\n" \
"	uint64 s0;\n" \
"	uint64 s1;\n" \
"} uint128;\n" \
"\n" \
"inline int ilog2_64(const uint64 x) { return 63 - clz(x); }\n" \
"\n" \
"inline int ilog2_96(const uint96 x) { return (x.s1 == 0) ? (63 - clz(x.s0)) : (95 - clz(x.s1)); }\n" \
"\n" \
"inline bool bittest_32(const uint32 x, const int b) { return ((x & ((uint32)(1) << b)) != 0); }\n" \
"\n" \
"inline bool bittest_64(const uint64 x, const int b) { return ((x & ((uint64)(1) << b)) != 0); }\n" \
"\n" \
"inline bool bittest_96(const uint96 x, const int b) { return (b < 64) ? bittest_64(x.s0, b) : bittest_32(x.s1, b - 64); }\n" \
"\n" \
"inline uint128 uint128_add(const uint128 x, const uint128 y)\n" \
"{\n" \
"	const uint64 s0 = x.s0 + y.s0;\n" \
"	const uint32 c = (s0 < y.s0) ? 1 : 0;\n" \
"	uint128 r; r.s0 = s0; r.s1 = x.s1 + y.s1 + c;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint128 uint128_add_64(const uint128 x, const uint64 y)\n" \
"{\n" \
"	const uint64 s0 = x.s0 + y;\n" \
"	const uint32 c = (s0 < y) ? 1 : 0;\n" \
"	uint128 r; r.s0 = s0; r.s1 = x.s1 + c;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint128 uint128_sub_64(const uint128 x, const uint64 y)\n" \
"{\n" \
"	const uint32 c = (x.s0 < y) ? 1 : 0;\n" \
"	uint128 r; r.s0 = x.s0 - y; r.s1 = x.s1 - c;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint128 uint128_sub_mod(const uint64 x, const uint128 y, const uint128 m)\n" \
"{\n" \
"	const uint32 c = (x < y.s0) ? 1 : 0;\n" \
"	uint128 r; r.s0 = x - y.s0; r.s1 = -y.s1 - c;\n" \
"	if (r.s1 != 0) r = uint128_add(r, m);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint128 uint128_mul_64_32(const uint64 x, const uint32 y) { uint128 r; r.s0 = x * y; r.s1 = mul_hi(x, (uint64)(y)); return r; }\n" \
"\n" \
"inline uint128 uint128_mul_64_64(const uint64 x, const uint64 y) { uint128 r; r.s0 = x * y; r.s1 = mul_hi(x, y); return r; }\n" \
"\n" \
"inline uint128 uint128_shift(const uint128 x, const int s)\n" \
"{\n" \
"	uint128 r;\n" \
"	r.s0 = (x.s0 >> s) | (x.s1 << (64 - s));\n" \
"	r.s1 = x.s1 >> s;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint64 uint128_64_shift(const uint128 x, const int s) { return (x.s0 >> s) | (x.s1 << (64 - s)); }\n" \
"\n" \
"inline uint32 kn_mod(const uint64 k, const int n, const uint32 m)\n" \
"{\n" \
"	const uint32 f = (1u << n) % m;\n" \
"	return ((uint32)(k % m) * f + 1) % m;\n" \
"}\n" \
"\n" \
"// Modular SIMD arithmetic in Mathemagix, Joris van der Hoeven, Grégoire Lecerf, Guillaume Quintin, Chapter 2.2.\n" \
"// n = 63, t = n + 1 = 64\n" \
"\n" \
"inline int divmod61_exp(const uint64 p) { return ilog2_64(p) + 1; }	// 2^{r-1} <= p < 2^r\n" \
"\n" \
"inline uint64 divmod61_invert(const uint64 p, const int r)\n" \
"{\n" \
"	// p_inv = 2^{s + 64} / p: Newton's method\n" \
"\n" \
"	// 2^{r - 2 + n + 1 - r} = 2^62 <= p_inv < 2^{r - 2 + n + 1 - (r - 1)} = 2^63\n" \
"	const int s = r - 2;\n" \
"\n" \
"	uint64 p_inv = 0, np_inv = (uint64)((~(uint32)(0)) / (uint32)(p >> 32)) << s;\n" \
"	while (p_inv < np_inv)\n" \
"	{\n" \
"		p_inv = np_inv;\n" \
"		np_inv = 2 * p_inv - uint128_64_shift(uint128_mul_64_64(mul_hi(p_inv, p_inv), p), s);\n" \
"	}\n" \
"\n" \
"	uint128 t = uint128_mul_64_64(np_inv, p);\n" \
"	const uint64 th_hi = (uint64)(1) << s;\n" \
"	while (t.s1 >= th_hi)\n" \
"	{\n" \
"		t = uint128_sub_64(t, p);\n" \
"		--np_inv;\n" \
"	}\n" \
"\n" \
"	return np_inv;\n" \
"}\n" \
"\n" \
"inline void divmod61(const uint128 a, const uint64 p, const uint64 p_inv, const int r, uint64 * const p_quot, uint64 * const p_rem)\n" \
"{\n" \
"	const int s = r - 2;	// alpha = 2^{n - 2} = 2^61, a < alpha * p => a < p^2 if p < 2^61\n" \
"\n" \
"	// a < 2^{2r} => b < 2^{r + 2} < 2^64 if p < 2^62\n" \
"	const uint64 b = uint128_64_shift(a, s);\n" \
"	uint64 quot = mul_hi(b, p_inv);\n" \
"	uint64 rem = a.s0 - quot * p;\n" \
"	// The number of iteration is 1 (see Table 3)\n" \
"	if (rem >= p) { rem -= p; ++quot; }\n" \
"	*p_quot = quot; *p_rem = rem;\n" \
"}\n" \
"\n" \
"typedef struct _Mod_k\n" \
"{\n" \
"	uint64 rem;		// 0 <= rem < k\n" \
"	uint32 quot;	// 0 <= quot <= 2^n\n" \
"} Mod_k;\n" \
"\n" \
"inline Mod_k norm(const Mod_k a, const uint64 k, const int n)\n" \
"{\n" \
"	const uint32 N = (uint32)(1) << n;\n" \
"	// k*2^n + k*2^n = k*2^{n + 1} => k*(2^n - 1) + k - 1\n" \
"	const bool is2N = (a.quot == 2*N);\n" \
"	Mod_k b; b.quot = is2N ? (N - 1) : a.quot; b.rem = is2N ? (k - 1) : a.rem;\n" \
"	// k*2^n = -1 (mod p)\n" \
"	if ((b.quot >= N) && (b.rem > 0)) { b.quot -= N; b.rem -= 1; }\n" \
"	return b;\n" \
"}\n" \
"\n" \
"inline Mod_k add_mod(const Mod_k a, const Mod_k b, const uint64 k, const int n)\n" \
"{\n" \
"	Mod_k c; c.rem = a.rem + b.rem; c.quot = a.quot + b.quot;\n" \
"	if (c.rem >= k) { c.rem -= k; c.quot += 1; }\n" \
"	return norm(c, k, n);\n" \
"}\n" \
"\n" \
"inline Mod_k mul_mod(const Mod_k a, const Mod_k b, const uint64 k, const int n, const uint64 k_inv, const int k_r)\n" \
"{\n" \
"	const uint32 N = (uint32)(1) << n;\n" \
"\n" \
"	// 0 <= ab_h < k * 2^{2n}, 0 <= ab_l < k^2\n" \
"	const uint128 ab_l = uint128_mul_64_64(a.rem, b.rem);\n" \
"	const uint128 ab_x = uint128_add(uint128_mul_64_32(a.rem, b.quot), uint128_mul_64_32(b.rem, a.quot));\n" \
"	const uint128 ab_h = uint128_add(uint128_mul_64_64(a.quot * (uint64)(b.quot), k), ab_x); \n" \
"\n" \
"	// 0 <= ab_quot < k * 2^{2n}, 0 <= ab_rem < k\n" \
"	uint64 ab_l_quot, ab_rem; divmod61(ab_l, k, k_inv, k_r, &ab_l_quot, &ab_rem);\n" \
"	const uint128 ab_quot = uint128_add_64(ab_h, ab_l_quot);\n" \
"\n" \
"	// (ab_quot >> n) * k * 2^n = -(ab_quot >> n) (mod p), (ab_quot >> n) < k * 2^n\n" \
"	uint128 p; p.s0 = (k << n) | 1; p.s1 = k >> (64 - n);\n" \
"	const uint128 r_l = uint128_sub_mod(ab_rem,	uint128_shift(ab_quot, n), p);\n" \
"	const uint32 r_h = (uint32)(ab_quot.s0) & (N - 1);\n" \
"\n" \
"	uint64 r_quot, r_rem; divmod61(r_l, k, k_inv, k_r, &r_quot, &r_rem);\n" \
"\n" \
"	Mod_k c; c.rem = r_rem; c.quot = r_quot + r_h;\n" \
"	return norm(c, k, n);\n" \
"}\n" \
"\n" \
"// a^e mod p, left-to-right algorithm\n" \
"inline Mod_k pow_mod(const Mod_k a, const uint64 e, const uint64 k, const int n, const uint64 k_inv, const int k_r)\n" \
"{\n" \
"	Mod_k r = a;\n" \
"	for (int b = ilog2_64(e) - 1; b >= 0; --b)\n" \
"	{\n" \
"		r = mul_mod(r, r, k, n, k_inv, k_r);\n" \
"		if (bittest_64(e, b)) r = mul_mod(r, a, k, n, k_inv, k_r);\n" \
"	}\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// 2^{(p - 1)/2} ?= +/-1 mod p\n" \
"inline bool prp(const uint64 k, const int n)\n" \
"{\n" \
"	const int k_r = divmod61_exp(k);\n" \
"	const uint64 k_inv = divmod61_invert(k, k_r);\n" \
"\n" \
"	// e = (p - 1)/2 = k*2^{n - 1}\n" \
"	uint96 e; e.s0 = k << (n - 1); e.s1 = (uint32)(k >> (64 - (n - 1)));\n" \
"	int b = ilog2_96(e) - 1;\n" \
"	Mod_k r; r.quot = 0; r.rem = bittest_96(e, b) ? 8 : 4; --b;\n" \
"	for (; b >= 64; --b)\n" \
"	{\n" \
"		r = mul_mod(r, r, k, n, k_inv, k_r);\n" \
"		if (bittest_32(e.s1, b - 64)) r = add_mod(r, r, k, n);\n" \
"	}\n" \
"	for (; b >= 0; --b)\n" \
"	{\n" \
"		r = mul_mod(r, r, k, n, k_inv, k_r);\n" \
"		if (bittest_64(e.s0, b)) r = add_mod(r, r, k, n);\n" \
"	}\n" \
"	const uint32 N = (uint32)(1) << n;\n" \
"	return (((r.quot == 0) && (r.rem == 1)) || ((r.quot == N) && (r.rem == 0)));\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void generate_primes(__global uint * restrict const prime_count, __global ulong * restrict const k_vector,\n" \
"	__global ulong2 * restrict const q_vector, __global ulong2 * restrict const one_vector, const ulong i)\n" \
"{\n" \
"	const uint64 k = (i << log2GlobalWorkSize) | get_global_id(0);\n" \
"\n" \
"	const int n = gfn_n + 1;\n" \
"	if (prp(k, n))\n" \
"	{\n" \
"		const uint prime_index = atomic_inc(prime_count);\n" \
"		k_vector[prime_index] = k;\n" \
"	}\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void init_factors(__global const uint * restrict const prime_count, __global const ulong * restrict const k_vector,\n" \
"	__global /*const*/ ulong2 * restrict const q_vector, __global const ulong2 * restrict const one_vector,\n" \
"	__global const char * restrict const _kro_vector,\n" \
"	__global ulong2 * restrict const c_vector, __global ulong2 * restrict const a2k_vector)\n" \
"{\n" \
"	const size_t i = get_global_id(0);\n" \
"	if (i >= *prime_count) return;\n" \
"\n" \
"	const uint64 k = k_vector[i];\n" \
"	const int n = gfn_n + 1;\n" \
"\n" \
"	// p = 1 (mod 4). If a is odd then (a/p) = (p/a) = ({p mod a}/a)\n" \
"\n" \
"	uint32 a = 3;\n" \
"	if (kn_mod(k, n, 3) != 2)\n" \
"	{\n" \
"		a += 2;\n" \
"		if (_kro_vector[256 * ((5 - 3) / 2) + kn_mod(k, n, 5)] >= 0)\n" \
"		{\n" \
"			a += 2;\n" \
"			if (_kro_vector[256 * ((7 - 3) / 2) + kn_mod(k, n, 7)] >= 0)\n" \
"			{\n" \
"				a += 4;\n" \
"				while (a < 256)\n" \
"				{\n" \
"					if (_kro_vector[256 * ((a - 3) / 2) + kn_mod(k, n, a)] < 0) break;\n" \
"					a += 2;\n" \
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
"	const int k_r = divmod61_exp(k);\n" \
"	const uint64 k_inv = divmod61_invert(k, k_r);\n" \
"	Mod_k c; c.quot = 0; c.rem = a;\n" \
"	c = pow_mod(c, k, k, n, k_inv, k_r);\n" \
"	const Mod_k a2k = mul_mod(c, c, k, n, k_inv, k_r);\n" \
"\n" \
"	c_vector[i] = (ulong2)(c.rem, c.quot);\n" \
"	a2k_vector[i] = (ulong2)(a2k.rem, a2k.quot);\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void check_factors(__global const uint * restrict const prime_count,\n" \
"	__global const ulong * restrict const k_vector, __global const ulong2 * restrict const q_vector,\n" \
"	__global ulong2 * restrict const c_vector, __global const ulong2 * restrict const a2k_vector,\n" \
"	__global uint * restrict const factor_count, __global ulong2 * restrict const factor)\n" \
"{\n" \
"	const size_t i = get_global_id(0);\n" \
"	if (i >= *prime_count) return;\n" \
"\n" \
"	const uint64 k = k_vector[i];\n" \
"	const int n = gfn_n + 1;\n" \
"\n" \
"	const ulong2 c_val = c_vector[i];\n" \
"	Mod_k c; c.rem = c_val.s0; c.quot = (uint32)(c_val.s1);\n" \
"	const ulong2 a2k_val = a2k_vector[i];\n" \
"	Mod_k a2k; a2k.rem = a2k_val.s0; a2k.quot = (uint32)(a2k_val.s1);\n" \
"\n" \
"	const int k_r = divmod61_exp(k);\n" \
"	const uint64 k_inv = divmod61_invert(k, k_r);\n" \
"\n" \
"	const uint32 N = (uint32)(1) << n;\n" \
"	for (size_t i = 0; i < factors_loop; ++i)\n" \
"	{\n" \
"		const uint32 parity = ((c.quot & (uint32)(k)) ^ (uint32)(c.rem)) & 1u;\n" \
"		if (parity != 0)\n" \
"		{\n" \
"			c.quot = N - c.quot;\n" \
"			const bool borrow = (c.rem > 1);\n" \
"			c.rem =  1 - c.rem;\n" \
"			if (borrow) { c.rem += k; c.quot -= 1; }\n" \
"		}\n" \
"\n" \
"		if ((c.quot == 0) && (c.rem <= 2000000000u))\n" \
"		{\n" \
"			uint128 p; p.s0 = (k << n) | 1; p.s1 = k >> (64 - n);\n" \
"			const uint factor_index = atomic_inc(factor_count);\n" \
"			factor[factor_index] = (ulong2)(p.s0, p.s1 | (c.rem << 32));\n" \
"		}\n" \
"\n" \
"		c = mul_mod(c, a2k, k, n, k_inv, k_r);		// c = a^{(2*i + 1).k}\n" \
"	}\n" \
"\n" \
"	c_vector[i] = (ulong2)(c.rem, c.quot);\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void clear_primes(__global uint * restrict const prime_count)\n" \
"{\n" \
"	*prime_count = 0;\n" \
"}\n" \
"";
