/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#ifdef __NV_CL_C_VERSION
	#define PTX_ASM	1
#endif

typedef uint	uint32;
typedef ulong	uint64;

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

inline int ilog2_64(const uint64 x) { return 63 - clz(x); }

inline int ilog2_96(const uint96 x) { return (x.s1 == 0) ? (63 - clz(x.s0)) : (95 - clz(x.s1)); }

inline bool bittest_32(const uint32 x, const int b) { return ((x & ((uint32)(1) << b)) != 0); }

inline bool bittest_64(const uint64 x, const int b) { return ((x & ((uint64)(1) << b)) != 0); }

inline bool bittest_96(const uint96 x, const int b) { return (b < 64) ? bittest_64(x.s0, b) : bittest_32(x.s1, b - 64); }

inline uint128 uint128_add(const uint128 x, const uint128 y)
{
	const uint64 s0 = x.s0 + y.s0;
	const uint32 c = (s0 < y.s0) ? 1 : 0;
	uint128 r; r.s0 = s0; r.s1 = x.s1 + y.s1 + c;
	return r;
}

inline uint128 uint128_add_64(const uint128 x, const uint64 y)
{
	const uint64 s0 = x.s0 + y;
	const uint32 c = (s0 < y) ? 1 : 0;
	uint128 r; r.s0 = s0; r.s1 = x.s1 + c;
	return r;
}

inline uint128 uint128_sub_64(const uint128 x, const uint64 y)
{
	const uint32 c = (x.s0 < y) ? 1 : 0;
	uint128 r; r.s0 = x.s0 - y; r.s1 = x.s1 - c;
	return r;
}

inline uint128 uint128_sub_mod(const uint64 x, const uint128 y, const uint128 m)
{
	const uint32 c = (x < y.s0) ? 1 : 0;
	uint128 r; r.s0 = x - y.s0; r.s1 = -y.s1 - c;
	if (r.s1 != 0) r = uint128_add(r, m);
	return r;
}

inline uint128 uint128_mul_64_32(const uint64 x, const uint32 y) { uint128 r; r.s0 = x * y; r.s1 = mul_hi(x, (uint64)(y)); return r; }

inline uint128 uint128_mul_64_64(const uint64 x, const uint64 y) { uint128 r; r.s0 = x * y; r.s1 = mul_hi(x, y); return r; }

inline uint128 uint128_shift(const uint128 x, const int s)
{
	uint128 r;
	r.s0 = (x.s0 >> s) | (x.s1 << (64 - s));
	r.s1 = x.s1 >> s;
	return r;
}

inline uint64 uint128_64_shift(const uint128 x, const int s) { return (x.s0 >> s) | (x.s1 << (64 - s)); }

inline uint32 kn_mod(const uint64 k, const int n, const uint32 m)
{
	const uint32 f = (1u << n) % m;
	return ((uint32)(k % m) * f + 1) % m;
}

// Modular SIMD arithmetic in Mathemagix, Joris van der Hoeven, Grégoire Lecerf, Guillaume Quintin, Chapter 2.2.
// n = 63, t = n + 1 = 64

inline int divmod61_exp(const uint64 p) { return ilog2_64(p) + 1; }	// 2^{r-1} <= p < 2^r

inline uint64 divmod61_invert(const uint64 p, const int r)
{
	// p_inv = 2^{s + 64} / p: Newton's method

	// 2^{r - 2 + n + 1 - r} = 2^62 <= p_inv < 2^{r - 2 + n + 1 - (r - 1)} = 2^63
	const int s = r - 2;

	uint64 p_inv = 0, np_inv = (uint64)((~(uint32)(0)) / ((uint32)(p >> 32) + 1)) << s;
	while (p_inv < np_inv)
	{
		p_inv = np_inv;
		np_inv = 2 * p_inv - uint128_64_shift(uint128_mul_64_64(mul_hi(p_inv, p_inv), p), s);
	}

	uint128 t = uint128_mul_64_64(np_inv, p);
	const uint64 th_hi = (uint64)(1) << s;
	while (t.s1 >= th_hi)
	{
		t = uint128_sub_64(t, p);
		--np_inv;
	}

	return np_inv;
}

inline void divmod61(const uint128 a, const uint64 p, const uint64 p_inv, const int r, uint64 * const p_quot, uint64 * const p_rem)
{
	const int s = r - 2;	// alpha = 2^{n - 2} = 2^61, a < alpha * p => a < p^2 if p < 2^61

	// a < 2^{2r} => b < 2^{r + 2} < 2^64 if p < 2^62
	const uint64 b = uint128_64_shift(a, s);
	uint64 quot = mul_hi(b, p_inv);
	uint64 rem = a.s0 - quot * p;
	// The number of iteration is 1 (see Table 3)
	if (rem >= p) { rem -= p; ++quot; }
	*p_quot = quot; *p_rem = rem;
}

typedef struct _Mod_k
{
	uint64 rem;		// 0 <= rem < k
	uint32 quot;	// 0 <= quot <= 2^n
} Mod_k;

inline Mod_k norm(const Mod_k a, const uint64 k, const int n)
{
	const uint32 N = (uint32)(1) << n;
	// k*2^n + k*2^n = k*2^{n + 1} => k*(2^n - 1) + k - 1
	const bool is2N = (a.quot == 2*N);
	Mod_k b; b.quot = is2N ? (N - 1) : a.quot; b.rem = is2N ? (k - 1) : a.rem;
	// k*2^n = -1 (mod p)
	if ((b.quot >= N) && (b.rem > 0)) { b.quot -= N; b.rem -= 1; }
	return b;
}

inline Mod_k add_mod(const Mod_k a, const Mod_k b, const uint64 k, const int n)
{
	Mod_k c; c.rem = a.rem + b.rem; c.quot = a.quot + b.quot;
	if (c.rem >= k) { c.rem -= k; c.quot += 1; }
	return norm(c, k, n);
}

inline Mod_k mul_mod(const Mod_k a, const Mod_k b, const uint64 k, const int n, const uint64 k_inv, const int k_r)
{
	const uint32 N = (uint32)(1) << n;

	// 0 <= ab_h < k * 2^{2n}, 0 <= ab_l < k^2
	const uint128 ab_l = uint128_mul_64_64(a.rem, b.rem);
	const uint128 ab_x = uint128_add(uint128_mul_64_32(a.rem, b.quot), uint128_mul_64_32(b.rem, a.quot));
	const uint128 ab_h = uint128_add(uint128_mul_64_64(a.quot * (uint64)(b.quot), k), ab_x); 

	// 0 <= ab_quot < k * 2^{2n}, 0 <= ab_rem < k
	uint64 ab_l_quot, ab_rem; divmod61(ab_l, k, k_inv, k_r, &ab_l_quot, &ab_rem);
	const uint128 ab_quot = uint128_add_64(ab_h, ab_l_quot);

	// (ab_quot >> n) * k * 2^n = -(ab_quot >> n) (mod p), (ab_quot >> n) < k * 2^n
	uint128 p; p.s0 = (k << n) | 1; p.s1 = k >> (64 - n);
	const uint128 r_l = uint128_sub_mod(ab_rem,	uint128_shift(ab_quot, n), p);
	const uint32 r_h = (uint32)(ab_quot.s0) & (N - 1);

	uint64 r_quot, r_rem; divmod61(r_l, k, k_inv, k_r, &r_quot, &r_rem);

	Mod_k c; c.rem = r_rem; c.quot = r_quot + r_h;
	return norm(c, k, n);
}

// a^e mod p, left-to-right algorithm
inline Mod_k pow_mod(const Mod_k a, const uint64 e, const uint64 k, const int n, const uint64 k_inv, const int k_r)
{
	Mod_k r = a;
	for (int b = ilog2_64(e) - 1; b >= 0; --b)
	{
		r = mul_mod(r, r, k, n, k_inv, k_r);
		if (bittest_64(e, b)) r = mul_mod(r, a, k, n, k_inv, k_r);
	}
	return r;
}

// 2^{(p - 1)/2} ?= +/-1 mod p
inline bool prp(const uint64 k, const int n, const uint64 k_inv, const int k_r)
{
	// e = (p - 1)/2 = k*2^{n - 1}
	uint96 e; e.s0 = k << (n - 1); e.s1 = (uint32)(k >> (64 - (n - 1)));
	int b = ilog2_96(e) - 1;
	Mod_k r; r.quot = 0; r.rem = bittest_96(e, b) ? 8 : 4; --b;
	for (; b >= 64; --b)
	{
		r = mul_mod(r, r, k, n, k_inv, k_r);
		if (bittest_32(e.s1, b - 64)) r = add_mod(r, r, k, n);
	}
	for (; b >= 0; --b)
	{
		r = mul_mod(r, r, k, n, k_inv, k_r);
		if (bittest_64(e.s0, b)) r = add_mod(r, r, k, n);
	}
	const uint32 N = (uint32)(1) << n;
	return (((r.quot == 0) && (r.rem == 1)) || ((r.quot == N) && (r.rem == 0)));
}

__kernel
void generate_primes(__global uint * restrict const prime_count, __global ulong * restrict const k_vector,
	__global ulong2 * restrict const k_ext_vector, const ulong i)
{
	const uint64 k = (i << log2GlobalWorkSize) | get_global_id(0);

	const int k_r = divmod61_exp(k);
	const uint64 k_inv = divmod61_invert(k, k_r);
	if (prp(k, g_n, k_inv, k_r))
	{
		const uint prime_index = atomic_inc(prime_count);
		k_vector[prime_index] = k;
		k_ext_vector[prime_index] = (ulong2)(k_inv, (ulong)(k_r));
	}
}

__kernel
void init_factors(__global const uint * restrict const prime_count, __global const ulong * restrict const k_vector,
	__global const ulong2 * restrict const k_ext_vector, __global const char * restrict const _kro_vector,
	__global ulong2 * restrict const c_vector, __global ulong2 * restrict const a2k_vector)
{
	const size_t i = get_global_id(0);
	if (i >= *prime_count) return;

	const uint64 k = k_vector[i];
	const ulong2 k_ext = k_ext_vector[i];
	const uint64 k_inv = k_ext.s0;
	const int k_r = (int)(k_ext.s1);

	// p = 1 (mod 4). If a is odd then (a/p) = (p/a) = ({p mod a}/a)

	uint32 a = 3;
	if (kn_mod(k, g_n, 3) != 2)
	{
		a += 2;
		if (_kro_vector[256 * ((5 - 3) / 2) + kn_mod(k, g_n, 5)] >= 0)
		{
			a += 2;
			if (_kro_vector[256 * ((7 - 3) / 2) + kn_mod(k, g_n, 7)] >= 0)
			{
				a += 4;
				while (a < 256)
				{
					if (_kro_vector[256 * ((a - 3) / 2) + kn_mod(k, g_n, a)] < 0) break;
					a += 2;
				}
				if (a >= 256)
				{
					c_vector[i] = (ulong2)(0, 0);
					return;
				}
			}
		}
	}

	Mod_k c; c.quot = 0; c.rem = a;
	c = pow_mod(c, k, k, g_n, k_inv, k_r);
	const Mod_k a2k = mul_mod(c, c, k, g_n, k_inv, k_r);

	c_vector[i] = (ulong2)(c.rem, c.quot);
	a2k_vector[i] = (ulong2)(a2k.rem, a2k.quot);
}

__kernel
void check_factors(__global const uint * restrict const prime_count,
	__global const ulong * restrict const k_vector, __global const ulong2 * restrict const k_ext_vector,
	__global ulong2 * restrict const c_vector, __global const ulong2 * restrict const a2k_vector,
	__global uint * restrict const factor_count, __global ulong2 * restrict const factor)
{
	const size_t i = get_global_id(0);
	if (i >= *prime_count) return;

	const uint64 k = k_vector[i];
	const ulong2 k_ext = k_ext_vector[i];
	const uint64 k_inv = k_ext.s0;
	const int k_r = (int)(k_ext.s1);

	const ulong2 c_val = c_vector[i];
	Mod_k c; c.rem = c_val.s0; c.quot = (uint32)(c_val.s1);
	const ulong2 a2k_val = a2k_vector[i];
	Mod_k a2k; a2k.rem = a2k_val.s0; a2k.quot = (uint32)(a2k_val.s1);

	const uint32 N = (uint32)(1) << g_n;
	for (size_t i = 0; i < factors_loop; ++i)
	{
		const uint32 parity = ((c.quot & (uint32)(k)) ^ (uint32)(c.rem)) & 1u;
		if (parity != 0)
		{
			c.quot = N - c.quot;
			const bool borrow = (c.rem > 1);
			c.rem =  1 - c.rem;
			if (borrow) { c.rem += k; c.quot -= 1; }
		}

		if ((c.quot == 0) && (c.rem <= 2000000000u))
		{
			uint128 p; p.s0 = (k << g_n) | 1; p.s1 = k >> (64 - g_n);
			const uint factor_index = atomic_inc(factor_count);
			factor[factor_index] = (ulong2)(p.s0, p.s1 | (c.rem << 32));
		}

		c = mul_mod(c, a2k, k, g_n, k_inv, k_r);		// c = a^{(2*i + 1).k}
	}

	c_vector[i] = (ulong2)(c.rem, c.quot);
}

__kernel
void clear_primes(__global uint * restrict const prime_count)
{
	*prime_count = 0;
}
