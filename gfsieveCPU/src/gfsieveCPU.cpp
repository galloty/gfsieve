/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <chrono>
#include <sys/stat.h>
#if defined(_WIN32)
#include <Windows.h>
#else
#include <signal.h>
#endif

inline int jacobi(const uint64_t x, const uint64_t y)
{
	uint64_t m = x, n = y;

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
		const uint64_t t = n; n = m; m = t;

		m %= n;	// (m/n) = (m mod n / n)
	}

	return 0;	// x and y are not coprime
}

// Modular operations, slow implementation
class Mod
{
private:
	const uint64_t _n;

private:
	static constexpr int ilog2(const uint64_t x) { return 63 - __builtin_clzll(x); }

public:
	Mod(const uint64_t n) : _n(n) { }

	uint64_t add(const uint64_t a, const uint64_t b) const
	{
		const uint64_t c = (a >= _n - b) ? _n : 0;
		return a + b - c;
	}

	uint64_t sub(const uint64_t a, const uint64_t b) const
	{
		const uint64_t c = (a < b) ? _n : 0;
		return a - b + c;
	}

	uint64_t mul(const uint64_t a, const uint64_t b) const
	{
		return uint64_t((a * __uint128_t(b)) % _n);
	}

	// a^e mod n, left-to-right algorithm
	uint64_t pow(const uint64_t a, const uint64_t e) const
	{
		uint64_t r = a;
		for (int b = ilog2(e) - 1; b >= 0; --b)
		{
			r = mul(r, r);
			if ((e & (uint64_t(1) << b)) != 0) r = mul(r, a);
		}
		return r;
	}

	bool spsp(const uint64_t p) const
	{
		// n - 1 = 2^k * r
		uint64_t r = _n - 1;
		int k = 0;
		for (; r % 2 == 0; r /= 2) ++k;

		uint64_t x = pow(p, r);
		if (x == 1) return true;

		// Compute x^(2^i) for 0 <= i < n.  If any are -1, n is a p-spsp.
		for (; k > 0; --k)
		{
			if (x == _n - 1) return true;
			x = mul(x, x);
		}

		return false;
	}

	bool isprime() const
	{
		if (_n < 2) return false;
		if (_n % 2 == 0) return (_n == 2);
		if (_n < 9) return true;

		// see https://oeis.org/A014233

		if (!spsp(2)) return false;
		if (_n < 2047ull) return true;

		if (!spsp(3)) return false;
		if (_n < 1373653ull) return true;

		if (!spsp(5)) return false;
		if (_n < 25326001ull) return true;

		if (!spsp(7)) return false;
		if (_n < 3215031751ull) return true;

		if (!spsp(11)) return false;
		if (_n < 2152302898747ull) return true;

		if (!spsp(13)) return false;
		if (_n < 3474749660383ull) return true;

		if (!spsp(17)) return false;
		if (_n < 341550071728321ull) return true;

		if (!spsp(19)) return false;
		// if (_n < 341550071728321ull) return true;

		if (!spsp(23)) return false;
		if (_n < 3825123056546413051ull) return true;

		if (!spsp(29)) return false;
		// if (_n < 3825123056546413051ull) return true;

		if (!spsp(31)) return false;
		// if (_n < 3825123056546413051ull) return true;

		if (!spsp(37)) return false;
		return true;	// 318665857834031151167461
	}
};

// Peter L. Montgomery, Modular multiplication without trial division, Math. Comp.44 (1985), 519–521.
class MpArith
{
private:
	const uint64_t _p, _q;
	const uint64_t _one;	// 2^64 mod p
	const uint64_t _r2;		// (2^64)^2 mod p

private:
	static constexpr int ilog2(const uint64_t x) { return 63 - __builtin_clzll(x); }
	static constexpr uint64_t mul_hi(const uint64_t lhs, const uint64_t rhs) { return uint64_t((lhs * __uint128_t(rhs)) >> 64); }

	// p * p_inv = 1 (mod 2^64) (Newton's method)
	static uint64_t invert(const uint64_t p)
	{
		uint64_t p_inv = 1, prev = 0;
		while (p_inv != prev) { prev = p_inv; p_inv *= 2 - p * p_inv; }
		return p_inv;
	}

	uint64_t REDC(const __uint128_t t) const
	{
		const uint64_t mp = mul_hi(uint64_t(t) * _q, _p), t_hi = uint64_t(t >> 64), r = t_hi - mp;
		return (t_hi < mp) ? r + _p : r;
	}

	uint64_t REDCshort(const uint64_t t) const
	{
		const uint64_t mp = mul_hi(t * _q, _p);
		return (mp != 0) ? _p - mp : 0;
	}

	// Montgomery form of 2^64 is (2^64)^2
	uint64_t two_pow_64() const
	{
		uint64_t t = add(_one, _one); t = add(t, t);	// 4
		t = add(t, t); t = add(t, t);					// 16
		for (size_t i = 0; i < 4; ++i) t = mul(t, t);	// 16^{2^4} = 2^64
		return t;
	}

public:
	MpArith(const uint64_t p) : _p(p), _q(invert(p)), _one((-p) % p), _r2(two_pow_64()) { }

	uint64_t toMp(const uint64_t n) const { return mul(n, _r2); }
	uint64_t toInt(const uint64_t r) const { return REDCshort(r); }

	uint64_t one() const { return _one; }
	uint64_t two() const { return add(_one, _one); }
	uint64_t three() const { return add(two(), _one); }

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

	// a^e mod p, left-to-right algorithm
	uint64_t pow(const uint64_t a, const uint64_t e) const
	{
		uint64_t r = a;
		for (int b = ilog2(e) - 1; b >= 0; --b)
		{
			r = mul(r, r);
			if ((e & (uint64_t(1) << b)) != 0) r = mul(r, a);
		}
		return r;
	}

	// 2^(p - 1) ?= 1 mod p
	bool prp() const
	{
		const uint64_t e = _p - 1;
		int b = ilog2(e) - 1;
		uint64_t r = add(_one, _one); r = add(r, r);
		if ((e & (uint64_t(1) << b)) != 0) r = add(r, r);
		for (--b; b >= 0; --b)
		{
			r = mul(r, r);
			if ((e & (uint64_t(1) << b)) != 0) r = add(r, r);
		}
		return (r == _one);
	}
};

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
	ss << "gfsieveCPU 25.1.0 " << sysver << ssc.str() << std::endl;
	ss << "Copyright (c) 2020, Yves Gallot" << std::endl;
	ss << "gfsieve is free source code, under the MIT license." << std::endl << std::endl;
	return ss.str();
}

static std::string usage()
{
	std::ostringstream ss;
	ss << "Usage: gfsieveCPU <n> <b_min> <b_max>" << std::endl;
	ss << "  n is exponent: b^{2^n} + 1, 8 <= n <= 24." << std::endl;
	ss << "  range is [b_min; b_max]." << std::endl;
	return ss.str();
}

inline std::string uint2string(const uint64_t n)
{
	double f = double(n); std::string prefix;
	if (f >= 1e15) { f *= 1e-15; prefix = "P"; }
	else if (f >= 1e12) { f *= 1e-12; prefix = "T"; }
	else if (f >= 1e9) { f *= 1e-9; prefix = "G"; }
	else if (f >= 1e6) { f *= 1e-6; prefix = "M"; }
	else return std::to_string(n);
	return std::to_string(f).substr(0, 5) + prefix;
}

class Sieve
{
private:
	const int _n;
	const uint64_t _b_min, _b_max;
	std::vector<bool> _bsieve;
	uint64_t _p_min;
	inline static volatile bool _quit = false;
	std::string _filename;
	std::chrono::high_resolution_clock::time_point _display_time, _record_time;

private:
	static void quit(int) { _quit = true; }

#if defined(_WIN32)
	static BOOL WINAPI HandlerRoutine(DWORD) { quit(1); return TRUE; }
#endif

public:
	Sieve(const int n, const uint64_t b_min, const uint64_t b_max)
		: _n(n), _b_min(b_min), _b_max(b_max), _bsieve((b_max - b_min) / 2 + 1, false)
	{
		_p_min = 0;
		std::stringstream ss; ss << "_" << n << "_" << b_min << "_" << b_max;
		_filename = ss.str();

#if defined(_WIN32)
		SetConsoleCtrlHandler(HandlerRoutine, TRUE);
#else
		signal(SIGTERM, quit);
		signal(SIGINT, quit);
#endif
	}

	virtual ~Sieve() {}

private:
	size_t get_size() const { return _bsieve.size(); }
	size_t get_count() const { size_t count = 0; for (bool b : _bsieve) if (!b) count++; return count; }

	std::string get_sieve_filename() const { return "sv" + _filename + ".dat"; }
	std::string get_cand_filename() const { return "cand" + _filename + ".txt"; }

	// Rosetta Code, CRC-32, C
	static uint32_t rc_crc32(const uint32_t crc32, const char * const buf, const size_t len)
	{
		static uint32_t table[256];
		static bool have_table = false;
	
		// This check is not thread safe; there is no mutex
		if (!have_table)
		{
			// Calculate CRC table
			for (size_t i = 0; i < 256; ++i)
			{
				uint32_t rem = static_cast<uint32_t>(i);  // remainder from polynomial division
				for (size_t j = 0; j < 8; ++j)
				{
					if (rem & 1)
					{
						rem >>= 1;
						rem ^= 0xedb88320;
					}
					else rem >>= 1;
				}
				table[i] = rem;
			}
			have_table = true;
		}

		uint32_t crc = ~crc32;
		for (size_t i = 0; i < len; ++i)
		{
			const uint8_t octet = static_cast<uint8_t>(buf[i]);  // Cast to unsigned octet
			crc = (crc >> 8) ^ table[(crc & 0xff) ^ octet];
		}
		return ~crc;
	}

	void info() const
	{
		static const double C[24] = { 0, 1.372813463, 2.678963880, 2.092794130, 3.671432123, 3.612924486, 3.942741295, 3.108964582,
			7.434805998, 7.489066280, 8.019343498, 7.224596905, 8.425349878, 8.467885720, 8.009684535, 5.802658835, 11.19571423,
			11.00430059, 13.00784637, 13.0724, 14.5167, 16.0846, 17.4099, 17.1286 };
		// #candidates = e^-gamma * C_n * (b_max - b_min) / log(p_max)
		const size_t size = 2 * get_size(), count = get_count();
		const double expected = 0.561459483567 * C[_n] * (_b_max - _b_min) / std::log(_p_min);
		std::cout << "p = " << uint2string(_p_min) << ", remaining " << count << "/" << size << " candidates (" << count * 100.0 / size << "%), expected: " << expected * 100.0 / size << "%." << std::endl;
	}

	bool read()
	{
		bool success = false;

		std::ifstream file(get_sieve_filename(), std::ios::binary);
		if (!file.good()) return false;

		file.read(reinterpret_cast<char *>(&_p_min), sizeof(_p_min));
		uint64_t sieve_size;
		file.read(reinterpret_cast<char *>(&sieve_size), sizeof(sieve_size));
		std::vector<uint32_t> sieve(sieve_size);
		file.read(reinterpret_cast<char *>(sieve.data()), static_cast<std::streamsize>(sieve_size * sizeof(uint32_t)));
		uint32_t crc32f = 0;
		file.read(reinterpret_cast<char *>(&crc32f), sizeof(crc32f));

		uint32_t crc32 = 0;
		crc32 = rc_crc32(crc32, reinterpret_cast<const char *>(&_p_min), sizeof(_p_min));
		crc32 = rc_crc32(crc32, reinterpret_cast<const char *>(&sieve_size), sizeof(sieve_size));
		crc32 = rc_crc32(crc32, reinterpret_cast<const char *>(sieve.data()), sieve_size * sizeof(uint32_t));
		crc32 = ~crc32 ^ 0xa23777ac;

		if (file && (crc32 == crc32f))
		{
			success = true;
			_bsieve.assign(_bsieve.size(), false);
			for (const uint32_t i : sieve) _bsieve[i] = true;
		}

		file.close();

		if (!success)
		{
			std::ostringstream ss; ss << "Error reading file '" << get_sieve_filename() << "'.";
			throw std::runtime_error(ss.str());
		}

		std::cout << "Resuming from a checkpoint." << std::endl;
		info();
		return true;
	}

	void write(const bool cand) const
	{
		struct stat s;
		const std::string sieve_filename = get_sieve_filename(), old_filename = sieve_filename + ".old";
		std::remove(old_filename.c_str());
		if ((stat(sieve_filename.c_str(), &s) == 0) && (std::rename(sieve_filename.c_str(), old_filename.c_str()) != 0))	// file exists and cannot rename it
		{
			std::ostringstream ss; ss << "Error writing file '" << sieve_filename << "'.";
			throw std::runtime_error(ss.str());
		}

		std::vector<uint32_t> sieve;
		for (size_t i = 0, size = _bsieve.size(); i < size; ++i) if (_bsieve[i]) sieve.push_back(uint32_t(i));
		const uint64_t sieve_size = sieve.size();

		uint32_t crc32 = 0;
		crc32 = rc_crc32(crc32, reinterpret_cast<const char *>(&_p_min), sizeof(_p_min));
		crc32 = rc_crc32(crc32, reinterpret_cast<const char *>(&sieve_size), sizeof(sieve_size));
		crc32 = rc_crc32(crc32, reinterpret_cast<const char *>(sieve.data()), sieve_size * sizeof(uint32_t));
		crc32 = ~crc32 ^ 0xa23777ac;

		std::ofstream file(sieve_filename, std::ios::binary);
		file.write(reinterpret_cast<const char *>(&_p_min), sizeof(_p_min));
		file.write(reinterpret_cast<const char *>(&sieve_size), sizeof(sieve_size));
		file.write(reinterpret_cast<const char *>(sieve.data()), static_cast<std::streamsize>(sieve.size() * sizeof(uint32_t)));
		file.write(reinterpret_cast<const char *>(&crc32), sizeof(crc32));
		file.close();

		if (cand)
		{
			std::ofstream file(get_cand_filename());
			for (size_t i = 0, size = get_size(); i < size; ++i) if (!_bsieve[i]) file << _b_min + 2 * i << std::endl;
			file.close();
		}

		info();
	}

	bool init()
	{
		_display_time = _record_time = std::chrono::high_resolution_clock::now();
		return !_quit;
	}

	bool monitor()
	{
		const bool quit = _quit;
		auto now = std::chrono::high_resolution_clock::now();
		if (quit || (std::chrono::duration<double>(now - _display_time).count() > 1))
		{
			_display_time = now;
			std::cout << uint2string(_p_min) << "\r";
			if (quit || (std::chrono::duration<double>(now - _record_time).count() > 5 * 60))
			{
				_record_time = now;
				write(quit);
			}
		}
		return !quit;
	}

	bool check_root(const uint64_t b, const uint64_t b_min, const uint64_t b_max, const uint64_t p, const int n)
	{
		for (uint64_t s = b; s <= b_max; s += p)
		{
			if (s >= b_min)
			{
				if (!_bsieve[(s - b_min) / 2])
				{
					const Mod mod(p);
					if (mod.pow(s, 1 << n) == p - 1)
					{
						_bsieve[(s - b_min) / 2] = true;
						// std::ofstream resFile("res.txt", std::ios::app);
						// if (resFile.is_open())
						// {
						// 	resFile << p << " " << s << std::endl;
						// 	resFile.close();
						// }
					}
					// May fail if p is not prime
					else if (mod.isprime())
					{
						std::ostringstream ss; ss << "Calculation error (check): p = " << p << ", b = " << s << ".";
						throw std::runtime_error(ss.str());
					}
					else return false;
				}
			}
		}
		return true;
	}

public:
	void check()
	{
		const int n = _n;
		const uint64_t b_min = _b_min, b_max = _b_max;

		std::cout << "n = " << n << ", b in [" << b_min << ", " << b_max << "]." << std::endl;

		if (!read()) _p_min = (1ul << (n + 1)) + 1;

		const uint64_t k_min = _p_min >> (n + 1), k_max = uint64_t(-1) >> (n + 1);

		if (!init()) return;

		for (uint64_t k = k_min; k <= k_max; ++k)
		{
			const uint64_t p = (k << (n + 1)) + 1;

			if ((p % 3 == 0) || (p % 5 == 0) || (p % 7 == 0) || (p % 11 == 0) || (p % 13 == 0) || (p % 17 == 0)
				|| (p % 19 == 0) || (p % 23 == 0) || (p % 29 == 0) || (p % 31 == 0)) continue;

			const MpArith mp(p);

			if (mp.prp())
			{
				uint64_t a = 3, ma = mp.three();
				while (jacobi(a, p) != -1) { ++a; ma = mp.add(ma, mp.one()); if (a == 256) { std::cout << p << std::endl; break; } }

				const uint64_t c = mp.pow(ma, k), b2 = mp.mul(c, c);
				uint64_t b = mp.toInt(c);

				for (uint64_t i = 0; i < (uint64_t(1) << (n - 1)); ++i)
				{
					const uint64_t beven = (b % 2 == 0) ? b : p - b;
					if (!check_root(beven, b_min, b_max, p, n)) break;
					b = mp.mul(b, b2);
				}

				_p_min = p;
				if (!monitor()) return;
			}
		}
	}
};

int main(int argc, char * argv[])
{
	std::cout << header();
	std::cout << std::fixed << std::setprecision(3);

	if (argc != 4)
	{
		std::cout << usage();
		return EXIT_SUCCESS;
	}

	int n = 0;
	uint64_t b_min = 0, b_max = 0;
	try
	{
		n = std::atoi(argv[1]);
		b_min = std::stoull(argv[2]) & ~1ull;
		b_max = std::stoull(argv[3]) & ~1ull;
	}
	catch (...)
	{
		std::cout << usage();
		return EXIT_SUCCESS;
	}

	if ((n < 8) || (n > 24)) { std::cerr << "n must be in [8, 24]." << std::endl; return EXIT_FAILURE; }
	if (b_min < 2) b_min = 2;
	if (b_max < b_min) b_max = b_min;

	Sieve sieve(n, b_min, b_max);

	try
	{
		sieve.check();
	}
	catch (const std::runtime_error & e)
	{
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
