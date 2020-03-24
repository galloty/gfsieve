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

inline uint96 uint96_set_ui(const uint32 n) { uint96 r; r.s0 = n; r.s1 = 0; return r; }

inline uint96 uint96_set(const uint64 s0, const uint32 s1) { uint96 r; r.s0 = s0; r.s1 = s1; return r; }

inline bool uint96_is_even(const uint96 x) { return (x.s0 % 2 == 0); }

inline bool uint96_is_equal_ui(const uint96 x, const uint32 n) { const bool r = ((x.s0 == n) && (x.s1 == 0)); return r; }

inline bool uint96_is_less_or_equal_ui(const uint96 x, const uint32 n) { const bool r = ((x.s0 <= n) && (x.s1 == 0)); return r; }

inline bool uint96_is_greater_or_equal(const uint96 x, const uint96 y)
{
	const bool r = (x.s1 > y.s1) || ((x.s1 == y.s1) && (x.s0 >= y.s0));
	return r;
}

inline int uint96_log2(const uint96 x)
{
	return (x.s1 == 0) ? (63 - clz(x.s0)) : (95 - clz(x.s1));
}

inline uint96 uint96_or_ui(const uint96 x, const uint32 n)
{
	uint96 r; r.s0 = x.s0 | n; r.s1 = x.s1;
	return r;
}

inline uint96 uint96_add_ui(const uint96 x, const uint32 n)
{
	const uint64 s0 = x.s0 + n;
	const uint32 c = (s0 < n) ? 1 : 0;
	uint96 r; r.s0 = s0; r.s1 = x.s1 + c;
	return r;
}

inline uint96 uint96_add(const uint96 x, const uint96 y)
{
	const uint64 s0 = x.s0 + y.s0;
	const uint32 c = (s0 < y.s0) ? 1 : 0;
	uint96 r; r.s0 = s0; r.s1 = x.s1 + y.s1 + c;
	return r;
}

inline uint96 uint96_sub_ui(const uint96 x, const uint32 n)
{
	const uint32 c = (x.s0 < n) ? 1 : 0;
	uint96 r; r.s0 = x.s0 - n; r.s1 = x.s1 - c;
	return r;
}

inline uint96 uint96_sub(const uint96 x, const uint96 y)
{
	const uint32 c = (x.s0 < y.s0) ? 1 : 0;
	uint96 r; r.s0 = x.s0 - y.s0; r.s1 = x.s1 - y.s1 - c;
	return r;
}

inline uint96 uint96_mul_2exp(const uint96 x, const int s)
{
	uint96 r;
	r.s0 = x.s0 << s;
	r.s1 = (x.s1 << s) | (uint32)(x.s0 >> (64 - s));
	return r;
}

inline uint96 uint96_div_2exp(const uint96 x, const int s)
{
	uint96 r;
	r.s0 = (x.s0 >> s) | ((uint64)(x.s1) << (64 - s));
	r.s1 = x.s1 >> s;
	return r;
}

inline uint96 uint96_mul_ui(const uint96 x, const uint32 n)
{
	const uint64 a0 = (uint32)(x.s0) * (uint64)(n);
	const uint64 a1 = (x.s0 >> 32) * n + (a0 >> 32);
	uint96 r;
	r.s0 = (uint32)(a0) | (a1 << 32);
	r.s1 = x.s1 * n + (a1 >> 32);
	return r;
}

inline uint96 uint96_mul(const uint96 x, const uint96 y)
{
	const uint64 a00l = x.s0 * y.s0, a00h = mul_hi(x.s0, y.s0);
	const uint32 a01 = (uint32)(x.s0) * y.s1, a10 = x.s1 * (uint32)(y.s0);
	const uint32 sa1 = (uint32)(a00h) + a01 + a10;
	uint96 r; r.s0 = a00l; r.s1 = sa1;
	return r;
}

inline uint128 uint128_add_ui(const uint128 x, const uint64 n)
{
	const uint64 s0 = x.s0 + n;
	const uint32 c = (s0 < n) ? 1 : 0;
	uint128 r; r.s0 = s0; r.s1 = x.s1 + c;
	return r;
}

inline uint128 uint128_add(const uint128 x, const uint128 y)
{
	const uint64 s0 = x.s0 + y.s0;
	const uint32 c = (s0 < y.s0) ? 1 : 0;
	uint128 r; r.s0 = s0; r.s1 = x.s1 + y.s1 + c;
	return r;
}

inline uint128 uint128_mul_ui(const uint64 n, const uint64 m)
{
	uint128 r; r.s0 = n * m; r.s1 = mul_hi(n, m);
	return r;
}

inline uint96 uint96_mul_hi(const uint96 x, const uint96 y)
{
	const uint128 a01 = uint128_mul_ui(x.s0, y.s1), a10 = uint128_mul_ui(x.s1, y.s0);
	const uint64 a00h = mul_hi(x.s0, y.s0), a11 = x.s1 * (uint64)(y.s1);

	uint128 sa1 = uint128_add_ui(uint128_add(a01, a10), a00h);
	sa1.s1 += a11;

	uint96 r; r.s0 = (sa1.s0 >> 32) | (sa1.s1 << 32); r.s1 = (uint32)(sa1.s1 >> 32);
	return r;
}

inline uint32 uint96_mod_ui(const uint96 x, const uint32 n)	// n < 2^8
{
	const uint32 a72 = x.s1 >> 8;
	const uint32 a48 = ((x.s1 & ((1 << 8) - 1)) << 16) | (uint32)(x.s0 >> 48);
	const uint32 a24 = (uint32)(x.s0 >> 24) & ((1 << 24) - 1);
	const uint32 a0 = (uint32)(x.s0) & ((1 << 24) - 1);

	uint32 r = a72 % n;
	r = ((r << 24) | a48) % n;
	r = ((r << 24) | a24) % n;
	r = ((r << 24) | a0) % n;
	return r;
}

inline uint96 uint96_dup_mod(const uint96 x, const uint96 p)
{
	uint96 r = uint96_mul_2exp(x, 1);
	if (uint96_is_greater_or_equal(r, p)) r = uint96_sub(r, p);
	return r;
}

inline double2 uint96_to_dd(const uint96 x)
{
	double2 r;
	r.s0 = ldexp((double)(((uint64)x.s1 << 16) | (x.s0 >> 48)), 48); r.s1 = 0;
	r = dd_add_d(r, (double)(x.s0 & (((uint64)1 << 48) - 1)));
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

inline uint96 uint96_barrett_inv(const uint96 p, const int p_shift)
{
	double2 d; d.s0 = ldexp(1.0, 96 + p_shift); d.s1 = 0;
	d = dd_floor(dd_div(d, uint96_to_dd(p)));
	return dd_to_uint96(d);
}

inline uint96 uint96_shoup_inv(const uint96 x, const uint96 p)
{
	double2 d; d.s0 = ldexp((double)(((uint64)x.s1 << 16) | (x.s0 >> 48)), 48 + 96); d.s1 = 0;
	d = dd_add_d(d, ldexp((double)(x.s0 & (((uint64)1 << 48) - 1)), 96));
	d = dd_floor(dd_div(d, uint96_to_dd(p)));
	return dd_to_uint96(d);
}

// Barrett's product: let n = 95, r = ceil(log2(p)), p_shift = r - 2 = ceil(log2(p)) - 1, t = n + 1 = 96,
// p_inv = floor(2^(s + t) / p). Then the number of iterations h = 1.
// We must have x^2 < alpha.p with alpha < 2^(n-2).enf then p < 2^(n-2) = 2^93.
inline uint96 uint96_square_mod(const uint96 x, const uint96 p, const uint96 p_inv, const int p_shift)
{
	const uint128 a00 = uint128_mul_ui(x.s0, x.s0);
	const uint128 a01 = uint128_mul_ui(x.s0, x.s1 + (uint64)(x.s1));
	const uint64 a11 = x.s1 * (uint64)(x.s1);

	uint128 sa1 = uint128_add_ui(a01, a00.s1);
	sa1.s1 += a11;

	const uint64 a0 = a00.s0, a1 = sa1.s0, a2 = sa1.s1;

	uint96 q_p;
	if (p_shift < 64)
	{
		q_p.s0 = (a0 >> p_shift) | (a1 << (64 - p_shift));
		q_p.s1 = (uint32)(a1 >> p_shift) | (uint32)(a2 << (64 - p_shift));
	}
	else
	{
		const int s = p_shift - 64;
		q_p.s0 = (a1 >> s) | (a2 << (64 - s));
		q_p.s1 = (uint32)(a2 >> s);
	}

	q_p = uint96_mul(uint96_mul_hi(q_p, p_inv), p);

	uint96 r = uint96_set(a0, (uint32)(a1));
	r = uint96_sub(r, q_p);
	if (uint96_is_greater_or_equal(r, p)) r = uint96_sub(r, p);
	return r;
}

// Shoup’s modular multiplication: p < 2^95 and w_inv = floor(w.2^96/p).
inline uint96 uint96_mul_mod(const uint96 x, const uint96 w, const uint96 p, const uint96 w_inv)
{
	const uint96 q = uint96_mul_hi(x, w_inv);
	uint96 r = uint96_sub(uint96_mul(x, w), uint96_mul(q, p));
	if (uint96_is_greater_or_equal(r, p)) r = uint96_sub(r, p);
	return r;
}

inline uint96 uint96_two_powm(const uint96 e, const uint96 p, const uint96 p_inv, const int p_shift)
{
	// x = 2^e mod p, left-to-right algorithm
	uint96 r = uint96_set_ui(1);
	if (e.s1 != 0)
	{
		for (int b = 31 - clz(e.s1); b >= 0; --b)
		{
			r = uint96_square_mod(r, p, p_inv, p_shift);
			if ((e.s1 & ((uint32)(1) << b)) != 0) r = uint96_dup_mod(r, p);
		}
	}
	for (int b = (e.s1 == 0) ? (63 - clz(e.s0)) : 63; b >= 0; --b)
	{
		r = uint96_square_mod(r, p, p_inv, p_shift);
		if ((e.s0 & ((uint64)(1) << b)) != 0) r = uint96_dup_mod(r, p);
	}
	return r;
}

inline uint96 uint96_powm(const uint96 a, const uint96 e, const uint96 p, const uint96 p_inv, const int p_shift, const uint96 a_inv)
{
	// x = a^e mod p, left-to-right algorithm
	uint96 r = uint96_set_ui(1);
	if (e.s1 != 0)
	{
		for (int b = 31 - clz(e.s1); b >= 0; --b)
		{
			r = uint96_square_mod(r, p, p_inv, p_shift);
			if ((e.s1 & ((uint32)(1) << b)) != 0) r = uint96_mul_mod(r, a, p, a_inv);
		}
	}
	for (int b = (e.s1 == 0) ? (63 - clz(e.s0)) : 63; b >= 0; --b)
	{
		r = uint96_square_mod(r, p, p_inv, p_shift);
		if ((e.s0 & ((uint64)(1) << b)) != 0) r = uint96_mul_mod(r, a, p, a_inv);
	}
	return r;
}

__kernel
void check_primes(__global uint * restrict const prime_count, __global ulong2 * restrict const prime_vector, const ulong i)
{
	const uint96 k = uint96_set((i << log2_prime_size) | get_global_id(0), (uint32)(i >> (64 - log2_prime_size)));

	const uint96 pm1 = uint96_mul_2exp(k, gfn_n + 1);
	const uint96 p = uint96_or_ui(pm1, 1);
	const int p_shift = uint96_log2(p) - 1;
	const uint96 p_inv = uint96_barrett_inv(p, p_shift);

	// 2-prp
	const uint96 r1 = uint96_two_powm(pm1, p, p_inv, p_shift);
	if (uint96_is_equal_ui(r1, 1))
	{
		const uint prime_index = atomic_inc(prime_count);
		prime_vector[prime_index] = (ulong2)(p.s0, (ulong)(p.s1));
	}
}

__kernel
void init_factors(__global const uint * restrict const prime_count, __global const ulong2 * restrict const prime_vector,
	__global ulong2 * restrict const a2k_vector, __global ulong2 * restrict const a2k_inv_vector,
	__global ulong2 * restrict const c_vector)
{
	const size_t i = get_global_id(0);
	if (i >= *prime_count) return;

	const ulong2 prm = prime_vector[i];
	const uint96 p = uint96_set(prm.s0, (uint32)(prm.s1));

	const uint96 k = uint96_div_2exp(p, gfn_n + 1);
	const int p_shift = uint96_log2(p) - 1;
	const uint96 p_inv = uint96_barrett_inv(p, p_shift);

	uint32 a = 3;
	for (; a < 256; a += 2)
	{
		const uint32 p_moda = uint96_mod_ui(p, a);
		if (jacobi(p_moda, a) == -1) break;
	}
	if (a >= 256) return;	// error?

	const uint96 a_inv = uint96_shoup_inv(uint96_set_ui(a), p);
	const uint96 c = uint96_powm(uint96_set_ui(a), k, p, p_inv, p_shift, a_inv);
	c_vector[i] = (ulong2)(c.s0, c.s1);

	const uint96 a2k = uint96_square_mod(c, p, p_inv, p_shift);
	a2k_vector[i] = (ulong2)(a2k.s0, a2k.s1);
	const uint96 a2k_inv = uint96_shoup_inv(a2k, p);
	a2k_inv_vector[i] = (ulong2)(a2k_inv.s0, a2k_inv.s1);
}

__kernel
void check_factors(__global const uint * restrict const prime_count, __global const ulong2 * restrict const prime_vector,
	__global const ulong2 * restrict const a2k_vector, __global const ulong2 * restrict const a2k_inv_vector,
	__global ulong2 * restrict const c_vector, __global uint * restrict const factor_count, __global ulong2 * restrict const factor)
{
	const size_t i = get_global_id(0);
	if (i >= *prime_count) return;

	const ulong2 prm = prime_vector[i];
	const uint96 p = uint96_set(prm.s0, (uint32)(prm.s1));
	const ulong2 a2k_val = a2k_vector[i];
	const uint96 a2k = uint96_set(a2k_val.s0, (uint32)(a2k_val.s1));
	const ulong2 a2k_inv_val = a2k_inv_vector[i];
	const uint96 a2k_inv = uint96_set(a2k_inv_val.s0, (uint32)(a2k_inv_val.s1));
	const ulong2 c_val = c_vector[i];
	uint96 c = uint96_set(c_val.s0, (uint32)(c_val.s1));

	for (uint32 i = 0; i < 1024; ++i)
	{
		// c is a^{(2*i + 1).k}

		const uint96 pmc = uint96_sub(p, c);
		const uint96 b = uint96_is_even(c) ? c : pmc;
		if (uint96_is_less_or_equal_ui(b, 2000000000))
		{
			const uint factor_index = atomic_inc(factor_count);
			factor[factor_index] = (ulong2)(p.s0, p.s1 | (b.s0 << 32));
		}

		c = uint96_mul_mod(c, a2k, p, a2k_inv);
	}

	c_vector[i] = (ulong2)(c.s0, c.s1);
}

__kernel
void clear_primes(__global uint * restrict const prime_count)
{
	*prime_count = 0;
}
