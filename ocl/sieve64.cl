/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#if !defined (__OPENCL_VERSION__)
	#define DOUBLE_EXTENSION	1
#elif defined (cl_khr_fp64)
	#define DOUBLE_EXTENSION	1
	#if __OPENCL_VERSION__ < 120
		#pragma OPENCL EXTENSION cl_khr_fp64: enable
	#endif
#elif defined (cl_amd_fp64)
	#define DOUBLE_EXTENSION	1
	#pragma OPENCL EXTENSION cl_amd_fp64: enable
#endif

#if !defined (DOUBLE_EXTENSION)
	#error "Double precision floating point not supported by OpenCL implementation."
#endif

#pragma OPENCL FP_CONTRACT OFF

#ifdef __NV_CL_C_VERSION
	#define PTX_ASM	1
#endif

typedef uint	uint32;
typedef ulong	uint64;

inline int jacobi(const uint32 x, const uint32 y)
{
	uint32 m = x, n = y;

	int k = 1;
	while (m != 0)
	{
		// (2/n) = (-1)^((n^2-1)/8)
		bool odd = false;
		while (m % 2 == 0) { m /= 2; odd = !odd; }
		if (odd && (n % 8 != 1) && (n % 8 != 7)) k = -k;

		if (m == 1) return k;	// (1/n) = 1

		// (m/n)(n/m) = -1 iif m == n == 3 (mod 4)
		if ((m % 4 == 3) && (n % 4 == 3)) k = -k;
		const uint32 t = n; n = m; m = t;

		m %= n;	// (m/n) = (m mod n / n)
	}

	return n;	// x and y are not coprime, return their gcd
}

inline double2 quick_sum(const double2 x)
{
	const double s = x.s0 + x.s1;
	const double e = x.s1 - (s - x.s0);
	return (double2)(s, e);
}

inline double2 sum(const double a, const double b)
{
	const double s = a + b;
	const double bb = s - a;
	const double e = (a - (s - bb)) + (b - bb);
	return (double2)(s, e);
}

inline double2 diff(const double a, const double b)
{
	const double d = a - b;
	const double bb = d - a;
	const double e = (a - (d - bb)) - (b + bb);
	return (double2)(d, e);
}

inline double2 prod(const double a, const double b)
{
	const double p = a * b;
	const double e = fma(a, b, -p);
	return (double2)(p, e);
}

inline double2 dd_add_d(const double2 x, const double d)
{
	double2 s = sum(x.s0, d);
	s.s1 += x.s1;
	return quick_sum(s);
}

inline double2 dd_sub(const double2 x, const double2 y)
{
	double2 s = diff(x.s0, y.s0);
	s.s1 += x.s1; s.s1 -= y.s1;
	return quick_sum(s);
}

inline double2 dd_sub_d(const double2 x, const double d)
{
	double2 s = diff(x.s0, d);
	s.s1 += x.s1;
	return quick_sum(s);
}

inline double2 dd_mul_d(const double2 x, const double d)
{
	double2 p = prod(x.s0, d);
	p.s1 += x.s1 * d;
	return quick_sum(p);
}

inline double2 dd_div(const double2 x, const double2 y)
{
	double2 q;
	q.s0 = x.s0 / y.s0;
	double2 z = dd_sub(x, dd_mul_d(y, q.s0));
	q.s1 = z.s0 / y.s0;
	z = dd_sub(z, dd_mul_d(y, q.s1));
	const double q3 = z.s0 / y.s0;
	return dd_add_d(quick_sum(q), q3);
}

inline double2 dd_floor(const double2 x)
{
	double2 r = (double2)(floor(x.s0), 0.0);
	if (r.s0 == x.s0)
	{
		r.s1 = floor(x.s1);
		r = quick_sum(r);
	}
	return r;
}

typedef struct _uint96
{
	uint64 s0;
	uint32 s1;
} uint96;

typedef struct _uint128
{
	uint64 s0;
	uint64 s1;
} uint128;

inline uint96 uint96_set_ui(const uint64 n) { uint96 r; r.s0 = n; r.s1 = 0; return r; }

inline uint96 uint96_set(const uint64 s0, const uint32 s1) { uint96 r; r.s0 = s0; r.s1 = s1; return r; }

inline bool uint96_is_negative(const uint96 x) { return ((int)(x.s1) < 0); }

inline int uint64_log2(const uint64 x) { return 63 - clz(x); }

inline uint96 uint96_add(const uint96 x, const uint96 y)
{
	uint96 r;
#ifdef PTX_ASM
	asm volatile ("add.cc.u64 %0, %1, %2;" : "=l" (r.s0) : "l" (x.s0), "l" (y.s0));
	asm volatile ("addc.u32 %0, %1, %2;" : "=r" (r.s1) : "r" (x.s1), "r" (y.s1));
#else
	const uint64 s0 = x.s0 + y.s0;
	const uint32 c = (s0 < y.s0) ? 1 : 0;
	r.s0 = s0; r.s1 = x.s1 + y.s1 + c;
#endif
	return r;
}

inline uint96 uint96_sub(const uint96 x, const uint96 y)
{
	uint96 r;
#ifdef PTX_ASM
	asm volatile ("sub.cc.u64 %0, %1, %2;" : "=l" (r.s0) : "l" (x.s0), "l" (y.s0));
	asm volatile ("subc.u32 %0, %1, %2;" : "=r" (r.s1) : "r" (x.s1), "r" (y.s1));
#else
	const uint32 c = (x.s0 < y.s0) ? 1 : 0;
	r.s0 = x.s0 - y.s0; r.s1 = x.s1 - y.s1 - c;
#endif
	return r;
}

inline uint96 uint96_mul(const uint96 x, const uint96 y)
{
	const uint32 a0 = (uint32)(x.s0), a1 = (uint32)(x.s0 >> 32), a2 = x.s1;
	const uint32 b0 = (uint32)(y.s0), b1 = (uint32)(y.s0 >> 32), b2 = y.s1;
	uint32 c0 = a0 * b0, c1 = mul_hi(a0, b0), c2 = a1 * b1;
	uint96 r;
#ifdef PTX_ASM
	asm volatile ("mad.lo.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a0), "r" (b2), "r" (c2));
	asm volatile ("mad.lo.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a2), "r" (b0), "r" (c2));
	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a0), "r" (b1), "r" (c1));
	asm volatile ("madc.hi.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a0), "r" (b1), "r" (c2));
	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a1), "r" (b0), "r" (c1));
	asm volatile ("madc.hi.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a1), "r" (b0), "r" (c2));
	r.s0 = upsample(c1, c0); r.s1 = c2;
#else
	c2 += a0 * b2 + a2 * b0;
	const uint64 c12 = c1 + a0 * (uint64)(b1) + a1 * (uint64)(b0);
	r.s0 = (c12 << 32) | c0; r.s1 = c2 + (uint32)(c12 >> 32);
#endif
	return r;
}

// x < 2^95, y < 2^96
inline uint96 uint96_mul_hi_ceil(const uint96 x, const uint96 y)
{
	uint96 r;
	const uint32 a0 = (uint32)(x.s0), a1 = (uint32)(x.s0 >> 32), a2 = x.s1;
	const uint32 b0 = (uint32)(y.s0), b1 = (uint32)(y.s0 >> 32), b2 = y.s1;
#ifdef PTX_ASM
	uint32 c1 = mul_hi(a0, b0), c2, c3, c4, c5;

	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a0), "r" (b1), "r" (c1));
	asm volatile ("madc.hi.u32 %0, %1, %2, 0;" : "=r" (c2) : "r" (a0), "r" (b1));

	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a0), "r" (b2), "r" (c2));
	asm volatile ("madc.hi.u32 %0, %1, %2, 1;" : "=r" (c3) : "r" (a0), "r" (b2));

	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a1), "r" (b0), "r" (c1));
	asm volatile ("madc.hi.cc.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a1), "r" (b0), "r" (c2));
	asm volatile ("madc.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c3) : "r" (a1), "r" (b2), "r" (c3));
	asm volatile ("madc.hi.u32 %0, %1, %2, 0;" : "=r" (c4) : "r" (a1), "r" (b2));

	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a1), "r" (b1), "r" (c2));
	asm volatile ("madc.hi.cc.u32 %0, %1, %2, %3;" : "=r" (c3) : "r" (a1), "r" (b1), "r" (c3));
	asm volatile ("addc.u32 %0, %1, 0;" : "=r" (c4) : "r" (c4));

	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a2), "r" (b0), "r" (c2));
	asm volatile ("madc.hi.cc.u32 %0, %1, %2, %3;" : "=r" (c3) : "r" (a2), "r" (b0), "r" (c3));
	asm volatile ("madc.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c4) : "r" (a2), "r" (b2), "r" (c4));
	asm volatile ("madc.hi.u32 %0, %1, %2, 0;" : "=r" (c5) : "r" (a2), "r" (b2));

	asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c3) : "r" (a2), "r" (b1), "r" (c3));
	asm volatile ("madc.hi.cc.u32 %0, %1, %2, %3;" : "=r" (c4) : "r" (a2), "r" (b1), "r" (c4));
	asm volatile ("addc.u32 %0, %1, 0;" : "=r" (c5) : "r" (c5));

	r.s0 = upsample(c4, c3); r.s1 = c5;
#else
	uint64 c1 = mul_hi(a0, b0), c2 = mul_hi(a0, b1), c3 = mad_hi(a0, b2, (uint32)(1)), c4 = mul_hi(a1, b2), c5 = mul_hi(a2, b2);

	c1 += a0 * b1; c1 += a1 * b0;
	c2 += mul_hi(a1, b0); c2 += a1 * b1; c2 += a0 * b2; c2 += a2 * b0; c2 += (c1 >> 32);
	c3 += mul_hi(a1, b1); c3 += a1 * b2; c3 += a2 * b1; c3 += mul_hi(a2, b0); c3 += (c2 >> 32);
	c4 += a2 * b2; c4 += mul_hi(a2, b1); c4 += (c3 >> 32);
	c5 += (c4 >> 32);

	r.s0 = ((uint64)(c4) << 32) | (uint32)c3; r.s1 = (uint32)c5;
#endif
	return r;
}

inline uint96 uint96_reduce_neg(const uint96 x, const uint96 p)
{
	return uint96_is_negative(x) ? uint96_add(x, p) : x;
}

inline uint64 uint64_dup_mod(const uint64 x, const uint64 p)
{
	const uint64 t = x + x;
	return (t >= p) ? t - p : t;
}

inline double2 uint64_to_dd(const uint64 x)
{
	double2 r;
	r.s0 = ldexp((double)(x >> 32), 32); r.s1 = 0;
	r = dd_add_d(r, (double)(x & (((uint64)1 << 32) - 1)));
	return r;
}

inline uint96 dd_to_uint96(const double2 x)
{
	double xh = floor(ldexp(x.s0, -48));
	const double2 d = dd_sub_d(x, ldexp(xh, 48));
	double xl = d.s0;
	if (xl < 0) { xh -= 1; xl += ldexp(1.0, 48); }
	const uint64 h = (uint64)(xh);
	const uint96 r = uint96_set((uint64)(xl) | (h << 48), (uint32)(h >> 16));
	return r;
}

inline uint96 uint96_barrett_inv(const uint64 p, const int p_shift)
{
	double2 d; d.s0 = ldexp(1.0, 96 + p_shift); d.s1 = 0;
	d = dd_floor(dd_div(d, uint64_to_dd(p)));
	return dd_to_uint96(d);
}

inline uint64 uint64_shoup_inv(const uint64 x, const uint64 p, const int nbits)
{
	double2 d; d.s0 = ldexp((double)(x >> 32), 32 + nbits); d.s1 = 0;
	d = dd_add_d(d, ldexp((double)(x & (((uint64)1 << 32) - 1)), nbits));
	d = dd_floor(dd_div(d, uint64_to_dd(p)));
	return dd_to_uint96(d).s0;
}

// Barrett's product: let n = 95, r = ceil(log2(p)), p_shift = r - 2 = ceil(log2(p)) - 1, t = n + 1 = 96,
// p_inv = floor(2^(s + t) / p). Then the number of iterations h = 1.
// We must have x^2 < alpha.p with alpha = 2^(n-2). If p <= 2^(n-2) = 2^93 then x^2 < p^2 <= alpha.p.
inline uint64 uint64_square_mod(const uint64 x, const uint64 p, const uint96 p_inv, const int p_shift)
{
#ifdef PTX_ASM
	const uint32 a0 = (uint32)(x), a1 = (uint32)(x >> 32);
	uint32 b0 = a0 * a0, b1, b2, b3;

	asm volatile ("mad.lo.cc.u32 %0, %1, %2, 0;" : "=r" (b1) : "r" (a0), "r" (a1));
	asm volatile ("madc.hi.u32 %0, %1, %2, 0;" : "=r" (b2) : "r" (a0), "r" (a1));

	asm volatile ("add.cc.u32 %0, %1, %2;" : "=r" (b1) : "r" (b1), "r" (b1));
	asm volatile ("addc.u32 %0, %1, %2;" : "=r" (b2) : "r" (b2), "r" (b2));

	asm volatile ("mad.hi.cc.u32 %0, %1, %2, %3;" : "=r" (b1) : "r" (a0), "r" (a0), "r" (b1));
	asm volatile ("madc.lo.cc.u32 %0, %1, %2, %3;" : "=r" (b2) : "r" (a1), "r" (a1), "r" (b2));
	asm volatile ("madc.hi.u32 %0, %1, %2, 0;" : "=r" (b3) : "r" (a1), "r" (a1));

	const uint64 x2_0 = upsample(b1, b0), x2_1 = upsample(b3, b2);
#else
	const uint64 x2_0 = x * x, x2_1 = mul_hi(x, x);
#endif

	uint96 q_p;
	q_p.s0 = (x2_0 >> p_shift) | (x2_1 << (64 - p_shift));
	q_p.s1 = (uint32)(x2_1 >> p_shift);

	const uint96 q = uint96_mul_hi_ceil(q_p, p_inv), p96 = uint96_set_ui(p);
	return uint96_reduce_neg(uint96_sub(uint96_set(x2_0, (uint32)(x2_1)), uint96_mul(q, p96)), p96).s0;
}

// Shoup’s modular multiplication: p < 2^63 and w_inv = floor(w.2^64/p).
inline uint64 uint64_mul_mod(const uint64 x, const uint64 w, const uint64 p, const uint64 w_inv)
{
	const uint64 q = mul_hi(x, w_inv);
	const uint64 t = x * w - q * p;
	return (t >= p) ? t - p : t;
}

inline uint64 uint64_two_powm(const uint64 e, const uint64 p, const uint96 p_inv, const int p_shift)
{
	// x = 2^e mod p, left-to-right algorithm
	uint64 r = 2;
	for (int b = uint64_log2(e) - 1; b >= 0; --b)
	{
		r = uint64_square_mod(r, p, p_inv, p_shift);
		if ((e & ((uint64)(1) << b)) != 0) r = uint64_dup_mod(r, p);
	}
	return r;
}

inline uint64 uint64_powm(const uint64 a, const uint64 e, const uint64 p, const uint96 p_inv, const int p_shift, const uint64 a_inv)
{
	// x = a^e mod p, left-to-right algorithm
	uint64 r = a;
	for (int b = uint64_log2(e) - 1; b >= 0; --b)
	{
		r = uint64_square_mod(r, p, p_inv, p_shift);
		if ((e & ((uint64)(1) << b)) != 0) r = uint64_mul_mod(r, a, p, a_inv);
	}
	return r;
}

__kernel
void check_primes(__global uint * restrict const prime_count, __global ulong2 * restrict const prime_vector, const ulong i)
{
	const uint64 k = (i << log2GlobalWorkSize) | get_global_id(0);

	const uint64 p = (k << (gfn_n + 1)) | 1;
	const int p_shift = uint64_log2(p) - 1;
	const uint96 p_inv = uint96_barrett_inv(p, p_shift);

	// 2-prp
	uint64 r = uint64_two_powm(k, p, p_inv, p_shift);
	for (size_t i = 0; i < gfn_n + 1; ++i) r = uint64_square_mod(r, p, p_inv, p_shift);

	if (r == 1)
	{
		const uint prime_index = atomic_inc(prime_count);
		prime_vector[prime_index] = (ulong2)(p, (ulong)(0));
	}
}

__kernel
void init_factors(__global const uint * restrict const prime_count, __global const ulong2 * restrict const prime_vector,
	__global ulong2 * restrict const a2k_vector, __global ulong2 * restrict const a2k_inv_vector,
	__global ulong2 * restrict const c_vector)
{
	const size_t i = get_global_id(0);
	if (i >= *prime_count) return;

	const uint64 p = prime_vector[i].s0;

	const uint64 k = p >> (gfn_n + 1);
	const int p_shift = uint64_log2(p) - 1;
	const uint96 p_inv = uint96_barrett_inv(p, p_shift);

	uint32 a;
	if ((p % 3) == 2) { a = 3; }
	else
	{
		const uint32 pmod5 = p % 5;
		if ((pmod5 == 2) || (pmod5 == 3)) { a = 5; }
		else
		{
			const uint32 pmod7 = p % 7;
			if ((pmod7 == 3) || (pmod7 == 5) || (pmod7 == 6)) { a = 7; }
			else
			{
				for (a = 11; a < 256; a += 2)
				{
					const uint32 pmoda = p % a;
					if (jacobi(pmoda, a) == -1) break;
				}
				if (a >= 256) return;	// error?
			}
		}
	}

	const uint64 a_inv = uint64_shoup_inv(a, p, 64);
	const uint64 c = uint64_powm(a, k, p, p_inv, p_shift, a_inv);
	c_vector[i] = (ulong2)(c, 0);

	const uint64 a2k = uint64_square_mod(c, p, p_inv, p_shift);
	a2k_vector[i] = (ulong2)(a2k, 0);
	const uint64 a2k_inv = uint64_shoup_inv(a2k, p, 64);
	a2k_inv_vector[i] = (ulong2)(a2k_inv, 0);
}

__kernel
void check_factors(__global const uint * restrict const prime_count, __global const ulong2 * restrict const prime_vector,
	__global const ulong2 * restrict const a2k_vector, __global const ulong2 * restrict const a2k_inv_vector,
	__global ulong2 * restrict const c_vector, __global uint * restrict const factor_count, __global ulong2 * restrict const factor)
{
	const size_t i = get_global_id(0);
	if (i >= *prime_count) return;

	const uint64 p = prime_vector[i].s0;
	const uint64 a2k = a2k_vector[i].s0;
	const uint64 a2k_inv = a2k_inv_vector[i].s0;
	uint64 c = c_vector[i].s0;

	for (size_t i = 0; i < factors_loop; ++i)
	{
		if (c % 2 != 0) c = p - c;
		if (c <= 2000000000)
		{
			const uint factor_index = atomic_inc(factor_count);
			factor[factor_index] = (ulong2)(p, c << 32);
		}
		c = uint64_mul_mod(c, a2k, p, a2k_inv);		// c = a^{(2*i + 1).k}
	}

	c_vector[i].s0 = c;
}

__kernel
void clear_primes(__global uint * restrict const prime_count)
{
	*prime_count = 0;
}
