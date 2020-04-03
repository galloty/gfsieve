/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

static const char * const src_ocl_sieve64 = \
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
"#ifdef __NV_CL_C_VERSION\n" \
"	#define PTX_ASM	1\n" \
"#endif\n" \
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
"inline uint96 uint96_set_ui(const uint64 n) { uint96 r; r.s0 = n; r.s1 = 0; return r; }\n" \
"\n" \
"inline uint96 uint96_set(const uint64 s0, const uint32 s1) { uint96 r; r.s0 = s0; r.s1 = s1; return r; }\n" \
"\n" \
"inline bool uint96_is_negative(const uint96 x) { return ((int)(x.s1) < 0); }\n" \
"\n" \
"inline int uint64_log2(const uint64 x) { return 63 - clz(x); }\n" \
"\n" \
"inline uint96 uint96_add(const uint96 x, const uint96 y)\n" \
"{\n" \
"	uint96 r;\n" \
"#ifdef PTX_ASM\n" \
"	asm volatile (\"add.cc.u64 %0, %1, %2;\" : \"=l\" (r.s0) : \"l\" (x.s0), \"l\" (y.s0));\n" \
"	asm volatile (\"addc.u32 %0, %1, %2;\" : \"=r\" (r.s1) : \"r\" (x.s1), \"r\" (y.s1));\n" \
"#else\n" \
"	const uint64 s0 = x.s0 + y.s0;\n" \
"	const uint32 c = (s0 < y.s0) ? 1 : 0;\n" \
"	r.s0 = s0; r.s1 = x.s1 + y.s1 + c;\n" \
"#endif\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_sub(const uint96 x, const uint96 y)\n" \
"{\n" \
"	uint96 r;\n" \
"#ifdef PTX_ASM\n" \
"	asm volatile (\"sub.cc.u64 %0, %1, %2;\" : \"=l\" (r.s0) : \"l\" (x.s0), \"l\" (y.s0));\n" \
"	asm volatile (\"subc.u32 %0, %1, %2;\" : \"=r\" (r.s1) : \"r\" (x.s1), \"r\" (y.s1));\n" \
"#else\n" \
"	const uint32 c = (x.s0 < y.s0) ? 1 : 0;\n" \
"	r.s0 = x.s0 - y.s0; r.s1 = x.s1 - y.s1 - c;\n" \
"#endif\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_mul(const uint96 x, const uint96 y)\n" \
"{\n" \
"	const uint32 a0 = (uint32)(x.s0), a1 = (uint32)(x.s0 >> 32), a2 = x.s1;\n" \
"	const uint32 b0 = (uint32)(y.s0), b1 = (uint32)(y.s0 >> 32), b2 = y.s1;\n" \
"	uint32 c0 = a0 * b0, c1 = mul_hi(a0, b0), c2 = a1 * b1;\n" \
"	uint96 r;\n" \
"#ifdef PTX_ASM\n" \
"	asm volatile (\"mad.lo.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a0), \"r\" (b2), \"r\" (c2));\n" \
"	asm volatile (\"mad.lo.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a2), \"r\" (b0), \"r\" (c2));\n" \
"	asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c1) : \"r\" (a0), \"r\" (b1), \"r\" (c1));\n" \
"	asm volatile (\"madc.hi.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a0), \"r\" (b1), \"r\" (c2));\n" \
"	asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c1) : \"r\" (a1), \"r\" (b0), \"r\" (c1));\n" \
"	asm volatile (\"madc.hi.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a1), \"r\" (b0), \"r\" (c2));\n" \
"	r.s0 = upsample(c1, c0); r.s1 = c2;\n" \
"#else\n" \
"	c2 += a0 * b2 + a2 * b0;\n" \
"	const uint64 c12 = c1 + a0 * (uint64)(b1) + a1 * (uint64)(b0);\n" \
"	r.s0 = (c12 << 32) | c0; r.s1 = c2 + (uint32)(c12 >> 32);\n" \
"#endif\n" \
"	return r;\n" \
"}\n" \
"\n" \
"// x < 2^95, y < 2^96\n" \
"inline uint96 uint96_mul_hi_ceil(const uint96 x, const uint96 y)\n" \
"{\n" \
"	uint96 r;\n" \
"	const uint32 a0 = (uint32)(x.s0), a1 = (uint32)(x.s0 >> 32), a2 = x.s1;\n" \
"	const uint32 b0 = (uint32)(y.s0), b1 = (uint32)(y.s0 >> 32), b2 = y.s1;\n" \
"#ifdef PTX_ASM\n" \
"	uint32 c1 = mul_hi(a0, b0), c2, c3, c4, c5;\n" \
"\n" \
"	asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c1) : \"r\" (a0), \"r\" (b1), \"r\" (c1));\n" \
"	asm volatile (\"madc.hi.u32 %0, %1, %2, 0;\" : \"=r\" (c2) : \"r\" (a0), \"r\" (b1));\n" \
"\n" \
"	asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a0), \"r\" (b2), \"r\" (c2));\n" \
"	asm volatile (\"madc.hi.u32 %0, %1, %2, 1;\" : \"=r\" (c3) : \"r\" (a0), \"r\" (b2));\n" \
"\n" \
"	asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c1) : \"r\" (a1), \"r\" (b0), \"r\" (c1));\n" \
"	asm volatile (\"madc.hi.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a1), \"r\" (b0), \"r\" (c2));\n" \
"	asm volatile (\"madc.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c3) : \"r\" (a1), \"r\" (b2), \"r\" (c3));\n" \
"	asm volatile (\"madc.hi.u32 %0, %1, %2, 0;\" : \"=r\" (c4) : \"r\" (a1), \"r\" (b2));\n" \
"\n" \
"	asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a1), \"r\" (b1), \"r\" (c2));\n" \
"	asm volatile (\"madc.hi.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c3) : \"r\" (a1), \"r\" (b1), \"r\" (c3));\n" \
"	asm volatile (\"addc.u32 %0, %1, 0;\" : \"=r\" (c4) : \"r\" (c4));\n" \
"\n" \
"	asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a2), \"r\" (b0), \"r\" (c2));\n" \
"	asm volatile (\"madc.hi.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c3) : \"r\" (a2), \"r\" (b0), \"r\" (c3));\n" \
"	asm volatile (\"madc.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c4) : \"r\" (a2), \"r\" (b2), \"r\" (c4));\n" \
"	asm volatile (\"madc.hi.u32 %0, %1, %2, 0;\" : \"=r\" (c5) : \"r\" (a2), \"r\" (b2));\n" \
"\n" \
"	asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c3) : \"r\" (a2), \"r\" (b1), \"r\" (c3));\n" \
"	asm volatile (\"madc.hi.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c4) : \"r\" (a2), \"r\" (b1), \"r\" (c4));\n" \
"	asm volatile (\"addc.u32 %0, %1, 0;\" : \"=r\" (c5) : \"r\" (c5));\n" \
"\n" \
"	r.s0 = upsample(c4, c3); r.s1 = c5;\n" \
"#else\n" \
"	uint64 c1 = mul_hi(a0, b0), c2 = mul_hi(a0, b1), c3 = mad_hi(a0, b2, (uint32)(1)), c4 = mul_hi(a1, b2), c5 = mul_hi(a2, b2);\n" \
"\n" \
"	c1 += a0 * b1; c1 += a1 * b0;\n" \
"	c2 += mul_hi(a1, b0); c2 += a1 * b1; c2 += a0 * b2; c2 += a2 * b0; c2 += (c1 >> 32);\n" \
"	c3 += mul_hi(a1, b1); c3 += a1 * b2; c3 += a2 * b1; c3 += mul_hi(a2, b0); c3 += (c2 >> 32);\n" \
"	c4 += a2 * b2; c4 += mul_hi(a2, b1); c4 += (c3 >> 32);\n" \
"	c5 += (c4 >> 32);\n" \
"\n" \
"	r.s0 = ((uint64)(c4) << 32) | (uint32)c3; r.s1 = (uint32)c5;\n" \
"#endif\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint96 uint96_reduce_neg(const uint96 x, const uint96 p)\n" \
"{\n" \
"	return uint96_is_negative(x) ? uint96_add(x, p) : x;\n" \
"}\n" \
"\n" \
"inline uint64 uint64_dup_mod(const uint64 x, const uint64 p)\n" \
"{\n" \
"	const uint64 t = x + x;\n" \
"	return (t >= p) ? t - p : t;\n" \
"}\n" \
"\n" \
"inline double2 uint64_to_dd(const uint64 x)\n" \
"{\n" \
"	double2 r;\n" \
"	r.s0 = ldexp((double)(x >> 32), 32); r.s1 = 0;\n" \
"	r = dd_add_d(r, (double)(x & (((uint64)1 << 32) - 1)));\n" \
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
"inline uint96 uint96_barrett_inv(const uint64 p, const int p_shift)\n" \
"{\n" \
"	double2 d; d.s0 = ldexp(1.0, 96 + p_shift); d.s1 = 0;\n" \
"	d = dd_floor(dd_div(d, uint64_to_dd(p)));\n" \
"	return dd_to_uint96(d);\n" \
"}\n" \
"\n" \
"inline uint64 uint64_shoup_inv(const uint64 x, const uint64 p, const int nbits)\n" \
"{\n" \
"	double2 d; d.s0 = ldexp((double)(x >> 32), 32 + nbits); d.s1 = 0;\n" \
"	d = dd_add_d(d, ldexp((double)(x & (((uint64)1 << 32) - 1)), nbits));\n" \
"	d = dd_floor(dd_div(d, uint64_to_dd(p)));\n" \
"	return dd_to_uint96(d).s0;\n" \
"}\n" \
"\n" \
"// Barrett's product: let n = 95, r = ceil(log2(p)), p_shift = r - 2 = ceil(log2(p)) - 1, t = n + 1 = 96,\n" \
"// p_inv = floor(2^(s + t) / p). Then the number of iterations h = 1.\n" \
"// We must have x^2 < alpha.p with alpha = 2^(n-2). If p <= 2^(n-2) = 2^93 then x^2 < p^2 <= alpha.p.\n" \
"inline uint64 uint64_square_mod(const uint64 x, const uint64 p, const uint96 p_inv, const int p_shift)\n" \
"{\n" \
"#ifdef PTX_ASM\n" \
"	const uint32 a0 = (uint32)(x), a1 = (uint32)(x >> 32);\n" \
"	uint32 b0 = a0 * a0, b1, b2, b3;\n" \
"\n" \
"	asm volatile (\"mad.lo.cc.u32 %0, %1, %2, 0;\" : \"=r\" (b1) : \"r\" (a0), \"r\" (a1));\n" \
"	asm volatile (\"madc.hi.u32 %0, %1, %2, 0;\" : \"=r\" (b2) : \"r\" (a0), \"r\" (a1));\n" \
"\n" \
"	asm volatile (\"add.cc.u32 %0, %1, %2;\" : \"=r\" (b1) : \"r\" (b1), \"r\" (b1));\n" \
"	asm volatile (\"addc.u32 %0, %1, %2;\" : \"=r\" (b2) : \"r\" (b2), \"r\" (b2));\n" \
"\n" \
"	asm volatile (\"mad.hi.cc.u32 %0, %1, %2, %3;\" : \"=r\" (b1) : \"r\" (a0), \"r\" (a0), \"r\" (b1));\n" \
"	asm volatile (\"madc.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (b2) : \"r\" (a1), \"r\" (a1), \"r\" (b2));\n" \
"	asm volatile (\"madc.hi.u32 %0, %1, %2, 0;\" : \"=r\" (b3) : \"r\" (a1), \"r\" (a1));\n" \
"\n" \
"	const uint64 x2_0 = upsample(b1, b0), x2_1 = upsample(b3, b2);\n" \
"#else\n" \
"	const uint64 x2_0 = x * x, x2_1 = mul_hi(x, x);\n" \
"#endif\n" \
"\n" \
"	uint96 q_p;\n" \
"	q_p.s0 = (x2_0 >> p_shift) | (x2_1 << (64 - p_shift));\n" \
"	q_p.s1 = (uint32)(x2_1 >> p_shift);\n" \
"\n" \
"	const uint96 q = uint96_mul_hi_ceil(q_p, p_inv), p96 = uint96_set_ui(p);\n" \
"	return uint96_reduce_neg(uint96_sub(uint96_set(x2_0, (uint32)(x2_1)), uint96_mul(q, p96)), p96).s0;\n" \
"}\n" \
"\n" \
"// Shoup’s modular multiplication: p < 2^63 and w_inv = floor(w.2^64/p).\n" \
"inline uint64 uint64_mul_mod(const uint64 x, const uint64 w, const uint64 p, const uint64 w_inv)\n" \
"{\n" \
"	const uint64 q = mul_hi(x, w_inv);\n" \
"	const uint64 t = x * w - q * p;\n" \
"	return (t >= p) ? t - p : t;\n" \
"}\n" \
"\n" \
"inline uint64 uint64_two_powm(const uint64 e, const uint64 p, const uint96 p_inv, const int p_shift)\n" \
"{\n" \
"	// x = 2^e mod p, left-to-right algorithm\n" \
"	uint64 r = 1;\n" \
"	for (int b = uint64_log2(e); b >= 0; --b)\n" \
"	{\n" \
"		r = uint64_square_mod(r, p, p_inv, p_shift);\n" \
"		if ((e & ((uint64)(1) << b)) != 0) r = uint64_dup_mod(r, p);\n" \
"	}\n" \
"	return r;\n" \
"}\n" \
"\n" \
"inline uint64 uint64_powm(const uint64 a, const uint64 e, const uint64 p, const uint96 p_inv, const int p_shift, const uint64 a_inv)\n" \
"{\n" \
"	// x = a^e mod p, left-to-right algorithm\n" \
"	uint64 r = 1;\n" \
"	for (int b = uint64_log2(e); b >= 0; --b)\n" \
"	{\n" \
"		r = uint64_square_mod(r, p, p_inv, p_shift);\n" \
"		if ((e & ((uint64)(1) << b)) != 0) r = uint64_mul_mod(r, a, p, a_inv);\n" \
"	}\n" \
"	return r;\n" \
"}\n" \
"\n" \
"__kernel\n" \
"void check_primes(__global uint * restrict const prime_count, __global ulong2 * restrict const prime_vector, const ulong i)\n" \
"{\n" \
"	const uint64 k = (i << log2GlobalWorkSize) | get_global_id(0);\n" \
"\n" \
"	const uint64 p = (k << (gfn_n + 1)) | 1;\n" \
"	const int p_shift = uint64_log2(p) - 1;\n" \
"	const uint96 p_inv = uint96_barrett_inv(p, p_shift);\n" \
"\n" \
"	// 2-prp\n" \
"	uint64 r = uint64_two_powm(k, p, p_inv, p_shift);\n" \
"	for (size_t i = 0; i < gfn_n + 1; ++i) r = uint64_square_mod(r, p, p_inv, p_shift);\n" \
"\n" \
"	if (r == 1)\n" \
"	{\n" \
"		const uint prime_index = atomic_inc(prime_count);\n" \
"		prime_vector[prime_index] = (ulong2)(p, (ulong)(0));\n" \
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
"	const uint64 p = prime_vector[i].s0;\n" \
"\n" \
"	const uint64 k = p >> (gfn_n + 1);\n" \
"	const int p_shift = uint64_log2(p) - 1;\n" \
"	const uint96 p_inv = uint96_barrett_inv(p, p_shift);\n" \
"\n" \
"	uint32 a;\n" \
"	if ((p % 3) == 2) { a = 3; }\n" \
"	else\n" \
"	{\n" \
"		const uint32 pmod5 = p % 5;\n" \
"		if ((pmod5 == 2) || (pmod5 == 3)) { a = 5; }\n" \
"		else\n" \
"		{\n" \
"			const uint32 pmod7 = p % 7;\n" \
"			if ((pmod7 == 3) || (pmod7 == 5) || (pmod7 == 6)) { a = 7; }\n" \
"			else\n" \
"			{\n" \
"				for (a = 11; a < 256; a += 2)\n" \
"				{\n" \
"					const uint32 pmoda = p % a;\n" \
"					if (jacobi(pmoda, a) == -1) break;\n" \
"				}\n" \
"				if (a >= 256) return;	// error?\n" \
"			}\n" \
"		}\n" \
"	}\n" \
"\n" \
"	const uint64 a_inv = uint64_shoup_inv(a, p, 64);\n" \
"	const uint64 c = uint64_powm(a, k, p, p_inv, p_shift, a_inv);\n" \
"	c_vector[i] = (ulong2)(c, 0);\n" \
"\n" \
"	const uint64 a2k = uint64_square_mod(c, p, p_inv, p_shift);\n" \
"	a2k_vector[i] = (ulong2)(a2k, 0);\n" \
"	const uint64 a2k_inv = uint64_shoup_inv(a2k, p, 64);\n" \
"	a2k_inv_vector[i] = (ulong2)(a2k_inv, 0);\n" \
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
"	const uint64 p = prime_vector[i].s0;\n" \
"	const uint64 a2k = a2k_vector[i].s0;\n" \
"	const uint64 a2k_inv = a2k_inv_vector[i].s0;\n" \
"	uint64 c = c_vector[i].s0;\n" \
"\n" \
"	for (size_t i = 0; i < factors_loop; ++i)\n" \
"	{\n" \
"		if (c % 2 != 0) c = p - c;\n" \
"		if (c <= 2000000000)\n" \
"		{\n" \
"			const uint factor_index = atomic_inc(factor_count);\n" \
"			factor[factor_index] = (ulong2)(p, c << 32);\n" \
"		}\n" \
"		c = uint64_mul_mod(c, a2k, p, a2k_inv);		// c = a^{(2*i + 1).k}\n" \
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
