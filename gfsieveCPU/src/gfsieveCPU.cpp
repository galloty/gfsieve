/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>

// Peter L. Montgomery, Modular multiplication without trial division, Math. Comp.44 (1985), 519–521.
class MpArith
{
private:
	const uint64_t _p, _q;
	const uint64_t _one;	// 2^64 mod p
	const uint64_t _r2;		// (2^64)^2 mod p

private:
	static int ilog2(const uint64_t x) { return 63 - __builtin_clzll(x); }

	// p * p_inv = 1 (mod 2^64) (Newton's method)
	static uint64_t invert(const uint64_t p)
	{
		uint64_t p_inv = 1, prev = 0;
		while (p_inv != prev) { prev = p_inv; p_inv *= 2 - p * p_inv; }
		return p_inv;
	}

	uint64_t REDC(const __uint128_t t) const
	{
		const uint64_t m = uint64_t(t) * _q;
		const int64_t r = int64_t((t >> 64) - uint64_t((m * __uint128_t(_p)) >> 64));
		return (r < 0) ? uint64_t(r + _p) : uint64_t(r);
	}

	uint64_t two_pow_64() const
	{
		uint64_t t = add(_one, _one); t = add(t, t);	// 4
		for (size_t i = 0; i < 5; ++i) t = mul(t, t);	// 4^{2^5} = 2^64
		return uint64_t(t);
	}

public:
	MpArith(const uint64_t p) : _p(p), _q(invert(p)), _one((-p) % p), _r2(two_pow_64()) { }

	uint64_t toMp(const uint64_t n) const { return mul(n, _r2); }
	uint64_t toInt(const uint64_t r) const { return REDC(r); }

	uint64_t add(const uint64_t a, const uint64_t b) const
	{
		const uint64_t c = (a >= _p - b) ? _p : 0;
		return a + b - c;
	}

	uint64_t sub(const uint64_t a, const uint64_t b) const
	{
		const uint64_t c = (a < b) ? _p : 0;
		return a - b + c;
	}

	uint64_t mul(const uint64_t a, const uint64_t b) const
	{
		return REDC(a * __uint128_t(b));
	}

	uint64_t pow(const uint64_t a, const uint64_t e) const
	{
		// x = a^e mod p, left-to-right algorithm
		uint64_t ra = toMp(a), r = ra;
		for (int b = ilog2(e) - 1; b >= 0; --b)
		{
			r = mul(r, r);
			if ((e & (uint64_t(1) << b)) != 0) r = mul(ra, r);
		}
		return toInt(r);
	}
};

inline int ilog2(const uint64_t x) { return 63 - __builtin_clzll(x); }

inline uint64_t powm(const uint32_t a, const uint64_t e, const uint64_t p)
{
	// x = a^e mod p, left-to-right algorithm
	uint64_t r = a;
	for (int b = ilog2(e) - 1; b >= 0; --b)
	{
		r = uint64_t((r * __uint128_t(r)) % p);
		if ((e & (uint64_t(1) << b)) != 0) r = uint64_t((a * __uint128_t(r)) % p);
	}
	return r;
}

static std::string header()
{
	const char * const sysver =
#if defined(_WIN64)
		"win64";
#elif defined(_WIN32)
		"win32";
#elif defined(__linux__)
#ifdef __x86_64
		"linux64";
#else
		"linux32";
#endif
#elif defined(__APPLE__)
		"macOS";
#else
		"unknown";
#endif

	std::ostringstream ssc;
#if defined(__GNUC__)
	ssc << " gcc-" << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
#elif defined(__clang__)
	ssc << " clang-" << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
#endif

	std::ostringstream ss;
	ss << "gfsieveCPU 0.3.0 " << sysver << ssc.str() << std::endl;
	ss << "Copyright (c) 2020, Yves Gallot" << std::endl;
	ss << "gfsieveCPU is free source code, under the MIT license." << std::endl << std::endl;
	return ss.str();
}

static std::string usage()
{
	std::ostringstream ss;
	ss << "Usage: gfsieveCPU <n> <b_min> <b_max>" << std::endl;
	ss << "         n is GFN exponent: b^{2^n} + 1, 10 <= n <= 24" << std::endl;
	ss << "         b_min is the start of the b range" << std::endl;
	ss << "         b_max is the end of the b range" << std::endl;
	return ss.str();
}

int main(int argc, char * argv[])
{
	std::cout << header();

	if (argc != 4)
	{
		std::cout << usage();
		return EXIT_FAILURE;
	}

	const uint32_t n = (argc > 1) ? std::min(std::max(std::atoi(argv[1]), 10), 24) : 22;
	const uint32_t b_min = (argc > 2) ? (std::atoi(argv[2]) / 2) * 2 : 100000800;
	const uint32_t b_max = (argc > 3) ?(std::atoi(argv[3]) / 2) * 2 : 100000900;

	std::stringstream ss; ss << "gfsieveCPU_" << n << "_" << b_min << "_" << b_max;
	const std::string res_filename = ss.str() + ".res";
	const std::string cand_filename = ss.str() + ".cand";

	std::vector<bool> sieve((b_max - b_min) / 2 + 1, false);

	const uint64_t p_max = 1000000000000000ull;	// 10^15
	uint64_t k_max = p_max >> (n + 1);
	while (((k_max << (n + 1)) + 1) < p_max) ++k_max;

	std::cout << "GFN-" << n << ": for k = 1 to " << k_max << ", for p = " << (uint64_t(1) << (n + 1)) + 1 << " to " << (k_max << (n + 1)) + 1 << std::endl;

	for (uint64_t k = 1; k <= k_max; ++k)
	{
		const uint64_t p = (k << (n + 1)) + 1;

		if (p % 3 == 0) continue;
		if (p % 5 == 0) continue;
		if (p % 7 == 0) continue;
		if (p % 11 == 0) continue;
		if (p % 13 == 0) continue;
		if (p % 17 == 0) continue;
		if (p % 19 == 0) continue;
		if (p % 23 == 0) continue;
		if (p % 29 == 0) continue;
		if (p % 31 == 0) continue;

		MpArith mp(p);

		if (mp.pow(1 << 16, (p - 1) / 16) == 1)	// 2-prp
		{
			const uint64_t t = powm(1 << 16, (p - 1) / 16, p); if (t != 1) std::cerr << "Error" << std::endl;
			for (uint32_t b = b_min, i = 0; b <= b_max; b += 2, ++i)
			{
				if (!sieve[i])
				{
					const uint64_t r = mp.pow(b, 1 << n);
					if (r == p - 1)
					{
						const uint64_t rp = powm(b, 1 << n, p); if (r != rp) std::cerr << "Error" << std::endl;
						sieve[i] = true;
						std::ofstream resFile(res_filename, std::ios::app);
						if (!resFile.is_open()) return EXIT_FAILURE;
						resFile << k << " * 2^" << n + 1 << " + 1 | " << b << "^" << (uint32_t(1) << n) << " + 1" << std::endl;
						resFile.flush();
						resFile.close();
						std::cout << k << " * 2^" << n + 1 << " + 1 | " << b << "^" << (uint32_t(1) << n) << " + 1" << std::endl;
						std::cout << int(100.0 * k / double(k_max)) <<"%\r";
					}
				}
			}
		}
	}

	size_t count = 0;
	for (bool b : sieve) if (b) count++;
	std::cout << "Removed " << count << "/" << sieve.size() << " candidates." << std::endl;

	std::ofstream candFile(cand_filename, std::ios::app);
	for (uint32_t b = b_min, i = 0; b <= b_max; b += 2, ++i) if (!sieve[i]) candFile << b << std::endl;
	candFile.close();

	return EXIT_SUCCESS;
}