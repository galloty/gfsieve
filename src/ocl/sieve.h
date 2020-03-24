/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

static const char * const src_ocl_sieve = \
"/*\n" \
"Copyright 2020, Yves Gallot\n" \
"\n" \
"gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.\n" \
"Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.\n" \
"*/\n" \
"\n" \
"#if !defined (__OPENCL_VERSION__)\n" \
"	#define DOUBLE_EXTENSION	1\n" \
"#elif defined (cl_khr_fp64)\n" \
"	#define DOUBLE_EXTENSION	1\n" \
"	#if __OPENCL_VERSION__ < 120\n" \
"		#pragma OPENCL EXTENSION cl_khr_fp64: enable\n" \
"	#endif\n" \
"#elif defined (cl_amd_fp64)\n" \
"	#define DOUBLE_EXTENSION	1\n" \
"	#pragma OPENCL EXTENSION cl_amd_fp64: enable\n" \
"#endif\n" \
"\n" \
"#if !defined (DOUBLE_EXTENSION)\n" \
"	#error \"Double precision floating point not supported by OpenCL implementation.\"\n" \
"#endif\n" \
"\n" \
"#pragma OPENCL FP_CONTRACT OFF\n" \
"\n" \
"typedef uint	uint32;\n" \
"typedef ulong	uint64;\n" \
"\n" \
"inline int jacobi(const uint32 x, const uint32 y)\n" \
"{\n" \
"	uint32 m = x, n = y;\n" \
"\n" \
"	int k = 1;\n" \
"	while (m != 0)\n" \
"	{\n" \
"		// (2/n) = (-1)^((n^2-1)/8)\n" \
"		bool odd = false;\n" \
"		while (m % 2 == 0) { m /= 2; odd = !odd; }\n" \
"		if (odd && (n % 8 != 1) && (n % 8 != 7)) k = -k;\n" \
"\n" \
"		if (m == 1) return k;	// (1/n) = 1\n" \
"\n" \
"		// (m/n)(n/m) = -1 iif m == n == 3 (mod 4)\n" \
"		if ((m % 4 == 3) && (n % 4 == 3)) k = -k;\n" \
"		const uint32 t = n; n = m; m = t;\n" \
"\n" \
"		m %= n;	// (m/n) = (m mod n / n)\n" \
"	}\n" \
"\n" \
"	return n;	// x and y are not coprime, return their gcd\n" \
"}\n" \
"\n" \
"inline double2 quick_sum(const double2 x)\n" \
"{\n" \
"	const double s = x.s0 + x.s1;\n" \
"	const double e = x.s1 - (s - x.s0);\n" \
"	return (double2)(s, e);\n" \
"}\n" \
"\n" \
"inline double2 sum(const double a, const double b)\n" \
"{\n" \
"	const double s = a + b;\n" \
"	const double bb = s - a;\n" \
"	const double e = (a - (s - bb)) + (b - bb);\n" \
"	return (double2)(s, e);\n" \
"}\n" \
"\n" \
"inline double2 diff(const double a, const double b)\n" \
"{\n" \
"	const double d = a - b;\n" \
"	const double bb = d - a;\n" \
"	const double e = (a - (d - bb)) - (b + bb);\n" \
"	return (double2)(d, e);\n" \
"}\n" \
"\n" \
"inline double2 prod(const double a, const double b)\n" \
"{\n" \
"	const double p = a * b;\n" \
"	const double e = fma(a, b, -p);\n" \
"	return (double2)(p, e);\n" \
"}\n" \
"\n" \
"inline double2 dd_add_d(const double2 x, const double d)\n" \
"{\n" \
"	double2 s = sum(x.s0, d);\n" \
"	s.s1 += x.s1;\n" \
"	return quick_sum(s);\n" \
"}\n" \
"\n" \
"inline double2 dd_sub(const double2 x, const double2 y)\n" \
"{\n" \
"	double2 s = diff(x.s0, y.s0);\n" \
"	s.s1 += x.s1; s.s1 -= y.s1;\n" \
"	return quick_sum(s);\n" \
"}\n" \
"\n" \
"inline double2 dd_sub_d(const double2 x, const double d)\n" \
"{\n" \
"	double2 s = diff(x.s0, d);\n" \
"	s.s1 += x.s1;\n" \
"	return quick_sum(s);\n" \
"}\n" \
"\n" \
"inline double2 dd_mul_d(const double2 x, const double d)\n" \
"{\n" \
"	double2 p = prod(x.s0, d);\n" \
"	p.s1 += x.s1 * d;\n" \
"	return quick_sum(p);\n" \
"}\n" \
"\n" \
"inline double2 dd_div(const double2 x, const double2 y)\n" \
"{\n" \
"	double2 q;\n" \
"	q.s0 = x.s0 / y.s0;\n" \
"	double2 z = dd_sub(x, dd_mul_d(y, q.s0));\n" \
"	q.s1 = z.s0 / y.s0;\n" \
"	z = dd_sub(z, dd_mul_d(y, q.s1));\n" \
"	const double q3 = z.s0 / y.s0;\n" \
"	return dd_add_d(quick_sum(q), q3);\n" \
"}\n" \
"\n" \
"inline double2 dd_floor(const double2 x)\n" \
"{\n" \
"	double2 r = (double2)(floor(x.s0), 0.0);\n" \
"	if (r.s0 == x.s0)\n" \
"	{\n" \
"		r.s1 = floor(x.s1);\n" \
"		r = quick_sum(r);\n" \
"	}\n" \
"	return r;\n" \
"}\n" \
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
"inline uint96 uint96_set_ui(const uint32 n) { uint96 r; r.s0 = n; r.s1 = 0; return r; }\n" \
"\n" \
"inline uint96 uint96_set(const uint64 s0, const uint32 s1) { uint96 r; r.s0 = s0; r.s1 = s1; return r; }\n" \
"\n" \
"inline bool uint96_is_even(const uint96 x) { return (x.s0 % 2 == 0); }\n" \
"\n" \
"inline bool uint96_is_equal_ui(const uint96 x, const uint32 n) { const bool r = ((x.s0 == n) && (x.s1 == 0)); return r; }\n" \
"\n" \
"inline bool uint96_is_less_or_equal_ui(const uint96 x, const uint32 n) { const bool r = ((x.s0 <= n) && (x.s1 == 0)); return r; }\n" \
"\n" \
"inline bool uint96_is_greater_or_equal(const uint96 x, const uint96 y)\n" \
"{\n" \
"	const bool r = (x.s1 > y.s1) || ((x.s1 == y.s1) && (x.s0 >= y.s0));\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline int uint96_log2(const uint96 x)\n" \
"{\n" \
"	return (x.s1 == 0) ? (63 - clz(x.s0)) : (95 - clz(x.s1));\n" \
"}\n" \
"\n" \
"inline uint96 uint96_or_ui(const uint96 x, const uint32 n)\n" \
"{\n" \
"	uint96 r; r.s0 = x.s0 | n; r.s1 = x.s1;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_add_ui(const uint96 x, const uint32 n)\n" \
"{\n" \
"	const uint64 s0 = x.s0 + n;\n" \
"	const uint32 c = (s0 < n) ? 1 : 0;\n" \
"	uint96 r; r.s0 = s0; r.s1 = x.s1 + c;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_add(const uint96 x, const uint96 y)\n" \
"{\n" \
"	const uint64 s0 = x.s0 + y.s0;\n" \
"	const uint32 c = (s0 < y.s0) ? 1 : 0;\n" \
"	uint96 r; r.s0 = s0; r.s1 = x.s1 + y.s1 + c;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_sub_ui(const uint96 x, const uint32 n)\n" \
"{\n" \
"	const uint32 c = (x.s0 < n) ? 1 : 0;\n" \
"	uint96 r; r.s0 = x.s0 - n; r.s1 = x.s1 - c;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_sub(const uint96 x, const uint96 y)\n" \
"{\n" \
"	const uint32 c = (x.s0 < y.s0) ? 1 : 0;\n" \
"	uint96 r; r.s0 = x.s0 - y.s0; r.s1 = x.s1 - y.s1 - c;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_mul_2exp(const uint96 x, const int s)\n" \
"{\n" \
"	uint96 r;\n" \
"	r.s0 = x.s0 << s;\n" \
"	r.s1 = (x.s1 << s) | (uint32)(x.s0 >> (64 - s));\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_div_2exp(const uint96 x, const int s)\n" \
"{\n" \
"	uint96 r;\n" \
"	r.s0 = (x.s0 >> s) | ((uint64)(x.s1) << (64 - s));\n" \
"	r.s1 = x.s1 >> s;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_mul_ui(const uint96 x, const uint32 n)\n" \
"{\n" \
"	const uint64 a0 = (uint32)(x.s0) * (uint64)(n);\n" \
"	const uint64 a1 = (x.s0 >> 32) * n + (a0 >> 32);\n" \
"	uint96 r;\n" \
"	r.s0 = (uint32)(a0) | (a1 << 32);\n" \
"	r.s1 = x.s1 * n + (a1 >> 32);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_mul(const uint96 x, const uint96 y)\n" \
"{\n" \
"	const uint64 a00l = x.s0 * y.s0, a00h = mul_hi(x.s0, y.s0);\n" \
"	const uint32 a01 = (uint32)(x.s0) * y.s1, a10 = x.s1 * (uint32)(y.s0);\n" \
"	const uint32 sa1 = (uint32)(a00h) + a01 + a10;\n" \
"	uint96 r; r.s0 = a00l; r.s1 = sa1;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint128 uint128_add_ui(const uint128 x, const uint64 n)\n" \
"{\n" \
"	const uint64 s0 = x.s0 + n;\n" \
"	const uint32 c = (s0 < n) ? 1 : 0;\n" \
"	uint128 r; r.s0 = s0; r.s1 = x.s1 + c;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint128 uint128_add(const uint128 x, const uint128 y)\n" \
"{\n" \
"	const uint64 s0 = x.s0 + y.s0;\n" \
"	const uint32 c = (s0 < y.s0) ? 1 : 0;\n" \
"	uint128 r; r.s0 = s0; r.s1 = x.s1 + y.s1 + c;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint128 uint128_mul_ui(const uint64 n, const uint64 m)\n" \
"{\n" \
"	uint128 r; r.s0 = n * m; r.s1 = mul_hi(n, m);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_mul_hi(const uint96 x, const uint96 y)\n" \
"{\n" \
"	const uint128 a01 = uint128_mul_ui(x.s0, y.s1), a10 = uint128_mul_ui(x.s1, y.s0);\n" \
"	const uint64 a00h = mul_hi(x.s0, y.s0), a11 = x.s1 * (uint64)(y.s1);\n" \
"\n" \
"	uint128 sa1 = uint128_add_ui(uint128_add(a01, a10), a00h);\n" \
"	sa1.s1 += a11;\n" \
"\n" \
"	uint96 r; r.s0 = (sa1.s0 >> 32) | (sa1.s1 << 32); r.s1 = (uint32)(sa1.s1 >> 32);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint32 uint96_mod_ui(const uint96 x, const uint32 n)	// n < 2^8\n" \
"{\n" \
"	const uint32 a72 = x.s1 >> 8;\n" \
"	const uint32 a48 = ((x.s1 & ((1 << 8) - 1)) << 16) | (uint32)(x.s0 >> 48);\n" \
"	const uint32 a24 = (uint32)(x.s0 >> 24) & ((1 << 24) - 1);\n" \
"	const uint32 a0 = (uint32)(x.s0) & ((1 << 24) - 1);\n" \
"\n" \
"	uint32 r = a72 % n;\n" \
"	r = ((r << 24) | a48) % n;\n" \
"	r = ((r << 24) | a24) % n;\n" \
"	r = ((r << 24) | a0) % n;\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_dup_mod(const uint96 x, const uint96 p)\n" \
"{\n" \
"	uint96 r = uint96_mul_2exp(x, 1);\n" \
"	if (uint96_is_greater_or_equal(r, p)) r = uint96_sub(r, p);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline double2 uint96_to_dd(const uint96 x)\n" \
"{\n" \
"	double2 r;\n" \
"	r.s0 = ldexp((double)(((uint64)x.s1 << 16) | (x.s0 >> 48)), 48); r.s1 = 0;\n" \
"	r = dd_add_d(r, (double)(x.s0 & (((uint64)1 << 48) - 1)));\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 dd_to_uint96(const double2 x)\n" \
"{\n" \
"	double xh = floor(ldexp(x.s0, -48));\n" \
"	const double2 d = dd_sub_d(x, ldexp(xh, 48));\n" \
"	double xl = d.s0;\n" \
"	if (xl < 0) { xh -= 1; xl += ldexp(1.0, 48); }\n" \
"	const uint64 h = (uint64)(xh);\n" \
"	const uint96 r = uint96_set((uint64)(xl) | (h << 48), (uint32)(h >> 16));\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_barrett_inv(const uint96 p, const int p_shift)\n" \
"{\n" \
"	double2 d; d.s0 = ldexp(1.0, 96 + p_shift); d.s1 = 0;\n" \
"	d = dd_floor(dd_div(d, uint96_to_dd(p)));\n" \
"	return dd_to_uint96(d);\n" \
"}\n" \
"\n" \
"inline uint96 uint96_shoup_inv(const uint96 x, const uint96 p)\n" \
"{\n" \
"	double2 d; d.s0 = ldexp((double)(((uint64)x.s1 << 16) | (x.s0 >> 48)), 48 + 96); d.s1 = 0;\n" \
"	d = dd_add_d(d, ldexp((double)(x.s0 & (((uint64)1 << 48) - 1)), 96));\n" \
"	d = dd_floor(dd_div(d, uint96_to_dd(p)));\n" \
"	return dd_to_uint96(d);\n" \
"}\n" \
"\n" \
"// Barrett's product: let n = 95, r = ceil(log2(p)), p_shift = r - 2 = ceil(log2(p)) - 1, t = n + 1 = 96,\n" \
"// p_inv = floor(2^(s + t) / p). Then the number of iterations h = 1.\n" \
"// We must have x^2 < alpha.p with alpha < 2^(n-2).enf then p < 2^(n-2) = 2^93.\n" \
"inline uint96 uint96_square_mod(const uint96 x, const uint96 p, const uint96 p_inv, const int p_shift)\n" \
"{\n" \
"	const uint128 a00 = uint128_mul_ui(x.s0, x.s0);\n" \
"	const uint128 a01 = uint128_mul_ui(x.s0, x.s1 + (uint64)(x.s1));\n" \
"	const uint64 a11 = x.s1 * (uint64)(x.s1);\n" \
"\n" \
"	uint128 sa1 = uint128_add_ui(a01, a00.s1);\n" \
"	sa1.s1 += a11;\n" \
"\n" \
"	const uint64 a0 = a00.s0, a1 = sa1.s0, a2 = sa1.s1;\n" \
"\n" \
"	uint96 q_p;\n" \
"	if (p_shift < 64)\n" \
"	{\n" \
"		q_p.s0 = (a0 >> p_shift) | (a1 << (64 - p_shift));\n" \
"		q_p.s1 = (uint32)(a1 >> p_shift) | (uint32)(a2 << (64 - p_shift));\n" \
"	}\n" \
"	else\n" \
"	{\n" \
"		const int s = p_shift - 64;\n" \
"		q_p.s0 = (a1 >> s) | (a2 << (64 - s));\n" \
"		q_p.s1 = (uint32)(a2 >> s);\n" \
"	}\n" \
"\n" \
"	q_p = uint96_mul(uint96_mul_hi(q_p, p_inv), p);\n" \
"\n" \
"	uint96 r = uint96_set(a0, (uint32)(a1));\n" \
"	r = uint96_sub(r, q_p);\n" \
"	if (uint96_is_greater_or_equal(r, p)) r = uint96_sub(r, p);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// Shoup’s modular multiplication: p < 2^95 and w_inv = floor(w.2^96/p).\n" \
"inline uint96 uint96_mul_mod(const uint96 x, const uint96 w, const uint96 p, const uint96 w_inv)\n" \
"{\n" \
"	const uint96 q = uint96_mul_hi(x, w_inv);\n" \
"	uint96 r = uint96_sub(uint96_mul(x, w), uint96_mul(q, p));\n" \
"	if (uint96_is_greater_or_equal(r, p)) r = uint96_sub(r, p);\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_two_powm(const uint96 e, const uint96 p, const uint96 p_inv, const int p_shift)\n" \
"{\n" \
"	// x = 2^e mod p, left-to-right algorithm\n" \
"	uint96 r = uint96_set_ui(1);\n" \
"	if (e.s1 != 0)\n" \
"	{\n" \
"		for (int b = 31 - clz(e.s1); b >= 0; --b)\n" \
"		{\n" \
"			r = uint96_square_mod(r, p, p_inv, p_shift);\n" \
"			if ((e.s1 & ((uint32)(1) << b)) != 0) r = uint96_dup_mod(r, p);\n" \
"		}\n" \
"	}\n" \
"	for (int b = (e.s1 == 0) ? (63 - clz(e.s0)) : 63; b >= 0; --b)\n" \
"	{\n" \
"		r = uint96_square_mod(r, p, p_inv, p_shift);\n" \
"		if ((e.s0 & ((uint64)(1) << b)) != 0) r = uint96_dup_mod(r, p);\n" \
"	}\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_powm(const uint96 a, const uint96 e, const uint96 p, const uint96 p_inv, const int p_shift, const uint96 a_inv)\n" \
"{\n" \
"	// x = a^e mod p, left-to-right algorithm\n" \
"	uint96 r = uint96_set_ui(1);\n" \
"	if (e.s1 != 0)\n" \
"	{\n" \
"		for (int b = 31 - clz(e.s1); b >= 0; --b)\n" \
"		{\n" \
"			r = uint96_square_mod(r, p, p_inv, p_shift);\n" \
"			if ((e.s1 & ((uint32)(1) << b)) != 0) r = uint96_mul_mod(r, a, p, a_inv);\n" \
"		}\n" \
"	}\n" \
"	for (int b = (e.s1 == 0) ? (63 - clz(e.s0)) : 63; b >= 0; --b)\n" \
"	{\n" \
"		r = uint96_square_mod(r, p, p_inv, p_shift);\n" \
"		if ((e.s0 & ((uint64)(1) << b)) != 0) r = uint96_mul_mod(r, a, p, a_inv);\n" \
"	}\n" \
"	return r;\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void check_primes(__global uint * restrict const prime_count, __global ulong2 * restrict const prime_vector, const ulong i)\n" \
"{\n" \
"	const uint96 k = uint96_set((i << log2_prime_size) | get_global_id(0), (uint32)(i >> (64 - log2_prime_size)));\n" \
"\n" \
"	const uint96 pm1 = uint96_mul_2exp(k, gfn_n + 1);\n" \
"	const uint96 p = uint96_or_ui(pm1, 1);\n" \
"	const int p_shift = uint96_log2(p) - 1;\n" \
"	const uint96 p_inv = uint96_barrett_inv(p, p_shift);\n" \
"\n" \
"	// 2-prp\n" \
"	const uint96 r1 = uint96_two_powm(pm1, p, p_inv, p_shift);\n" \
"	if (uint96_is_equal_ui(r1, 1))\n" \
"	{\n" \
"		const uint prime_index = atomic_inc(prime_count);\n" \
"		prime_vector[prime_index] = (ulong2)(p.s0, (ulong)(p.s1));\n" \
"	}\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void init_factors(__global const uint * restrict const prime_count, __global const ulong2 * restrict const prime_vector,\n" \
"	__global ulong2 * restrict const a2k_vector, __global ulong2 * restrict const a2k_inv_vector,\n" \
"	__global ulong2 * restrict const c_vector)\n" \
"{\n" \
"	const size_t i = get_global_id(0);\n" \
"	if (i >= *prime_count) return;\n" \
"\n" \
"	const ulong2 prm = prime_vector[i];\n" \
"	const uint96 p = uint96_set(prm.s0, (uint32)(prm.s1));\n" \
"\n" \
"	const uint96 k = uint96_div_2exp(p, gfn_n + 1);\n" \
"	const int p_shift = uint96_log2(p) - 1;\n" \
"	const uint96 p_inv = uint96_barrett_inv(p, p_shift);\n" \
"\n" \
"	uint32 a = 3;\n" \
"	for (; a < 256; a += 2)\n" \
"	{\n" \
"		const uint32 p_moda = uint96_mod_ui(p, a);\n" \
"		if (jacobi(p_moda, a) == -1) break;\n" \
"	}\n" \
"	if (a >= 256) return;	// error?\n" \
"\n" \
"	const uint96 a_inv = uint96_shoup_inv(uint96_set_ui(a), p);\n" \
"	const uint96 c = uint96_powm(uint96_set_ui(a), k, p, p_inv, p_shift, a_inv);\n" \
"	c_vector[i] = (ulong2)(c.s0, c.s1);\n" \
"\n" \
"	const uint96 a2k = uint96_square_mod(c, p, p_inv, p_shift);\n" \
"	a2k_vector[i] = (ulong2)(a2k.s0, a2k.s1);\n" \
"	const uint96 a2k_inv = uint96_shoup_inv(a2k, p);\n" \
"	a2k_inv_vector[i] = (ulong2)(a2k_inv.s0, a2k_inv.s1);\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void check_factors(__global const uint * restrict const prime_count, __global const ulong2 * restrict const prime_vector,\n" \
"	__global const ulong2 * restrict const a2k_vector, __global const ulong2 * restrict const a2k_inv_vector,\n" \
"	__global ulong2 * restrict const c_vector, __global uint * restrict const factor_count, __global ulong2 * restrict const factor)\n" \
"{\n" \
"	const size_t i = get_global_id(0);\n" \
"	if (i >= *prime_count) return;\n" \
"\n" \
"	const ulong2 prm = prime_vector[i];\n" \
"	const uint96 p = uint96_set(prm.s0, (uint32)(prm.s1));\n" \
"	const ulong2 a2k_val = a2k_vector[i];\n" \
"	const uint96 a2k = uint96_set(a2k_val.s0, (uint32)(a2k_val.s1));\n" \
"	const ulong2 a2k_inv_val = a2k_inv_vector[i];\n" \
"	const uint96 a2k_inv = uint96_set(a2k_inv_val.s0, (uint32)(a2k_inv_val.s1));\n" \
"	const ulong2 c_val = c_vector[i];\n" \
"	uint96 c = uint96_set(c_val.s0, (uint32)(c_val.s1));\n" \
"\n" \
"	for (uint32 i = 0; i < 1024; ++i)\n" \
"	{\n" \
"		// c is a^{(2*i + 1).k}\n" \
"\n" \
"		const uint96 pmc = uint96_sub(p, c);\n" \
"		const uint96 b = uint96_is_even(c) ? c : pmc;\n" \
"		if (uint96_is_less_or_equal_ui(b, 2000000000))\n" \
"		{\n" \
"			const uint factor_index = atomic_inc(factor_count);\n" \
"			factor[factor_index] = (ulong2)(p.s0, p.s1 | (b.s0 << 32));\n" \
"		}\n" \
"\n" \
"		c = uint96_mul_mod(c, a2k, p, a2k_inv);\n" \
"	}\n" \
"\n" \
"	c_vector[i] = (ulong2)(c.s0, c.s1);\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void clear_primes(__global uint * restrict const prime_count)\n" \
"{\n" \
"	*prime_count = 0;\n" \
"}\n" \
"";
