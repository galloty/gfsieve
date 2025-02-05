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
#include <queue>
#include <thread>
#include <mutex>
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
	uint64_t _p, _q;
	uint64_t _one, _r2;	// 2^64 mod p and (2^64)^2 mod p

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
	MpArith() : _p(0), _q(0), _one(0), _r2(0) {}
	MpArith(const uint64_t p) : _p(p), _q(invert(p)), _one((-p) % p), _r2(two_pow_64()) {}
	MpArith(const MpArith & rhs) : _p(rhs._p), _q(rhs._q), _one(rhs._one), _r2(rhs._r2) {}
	MpArith & operator=(const MpArith & rhs) { _p = rhs._p; _q = rhs._q; _one = rhs._one, _r2 = rhs._r2; return *this; }

	uint64_t p() const { return _p; }
	uint64_t q() const { return _q; }

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

using vu64 = uint64_t[4];

class MpArith4
{
private:
	vu64 _p, _q;

private:
	static constexpr uint64_t mul_hi(const uint64_t lhs, const uint64_t rhs) { return uint64_t((lhs * __uint128_t(rhs)) >> 64); }

public:
	MpArith4(const vu64 p, const vu64 q) { for (size_t i = 0; i < 4; ++i) { _p[i] = p[i]; _q[i] = q[i]; } }

	void mul(vu64 a, const vu64 b) const
	{
		for (size_t i = 0; i < 4; ++i)
		{
			const __uint128_t t = a[i] * __uint128_t(b[i]);
			const uint64_t mp = mul_hi(uint64_t(t) * _q[i], _p[i]), t_hi = uint64_t(t >> 64), r = t_hi - mp;
			a[i] = (t_hi < mp) ? r + _p[i] : r;
		}
	}
};

template<typename T>
class Fifo
{
private:
	static const size_t max_queue_size = 128;

	std::mutex _mutex;
	std::queue<T> _queue;
	bool _end = false;

public:
	void end() { std::lock_guard<std::mutex> guard(_mutex); _end = true; }

	void push(const T & val)
	{
		_mutex.lock();
		while (_queue.size() >= max_queue_size)
		{
			_mutex.unlock();
			std::this_thread::sleep_for(std::chrono::milliseconds(100));
			_mutex.lock();
		}
		_queue.push(val);
		_mutex.unlock();
	}

	bool pop(T & val)
	{
		_mutex.lock();
		while (_queue.empty())
		{
			if (_end)
			{
				_mutex.unlock();
				return false;
			}
			_mutex.unlock();
			std::cout << "Warning: empty FIFO." << std::endl;
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
			_mutex.lock();
		}
		val = _queue.front();
		_queue.pop();
		_mutex.unlock();
		return true;
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

class Sieve
{
private:
	const int _n;
	const uint64_t _b_min, _b_max;
	std::vector<bool> _bsieve;
	uint64_t _p_min;
	inline static volatile bool _quit = false;
	std::string _filename;
	std::chrono::high_resolution_clock::time_point _display_time, _record_time, _start_time;

	static const size_t karray_size = 1024;
	struct kArray { uint64_t k[karray_size]; };

	static const size_t parray_size = 1024;
	struct Mpbb2 { uint64_t p; uint64_t q; uint64_t b; uint64_t b2; };
	struct pArray { Mpbb2 mpbb2[parray_size]; };

	Fifo<kArray> _kqueue;	//	512KB
	Fifo<pArray> _pqueue;	//	2MB

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

	static std::string uint2string(const uint64_t n)
	{
		double f = double(n); std::string prefix;
		if (f >= 1e15) { f *= 1e-15; prefix = "P"; }
		else if (f >= 1e12) { f *= 1e-12; prefix = "T"; }
		else if (f >= 1e9) { f *= 1e-9; prefix = "G"; }
		else if (f >= 1e6) { f *= 1e-6; prefix = "M"; }
		else return std::to_string(n);
		return std::to_string(f).substr(0, 5) + prefix;
	}

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
		_display_time = _record_time = _start_time = std::chrono::high_resolution_clock::now();
		return !_quit;
	}

	bool monitor()
	{
		const bool quit = _quit;
		const auto now = std::chrono::high_resolution_clock::now();
		if (quit || (std::chrono::duration<double>(now - _display_time).count() > 1))
		{
			_display_time = now;
			std::cout << std::lrint(std::chrono::duration<double>(now - _start_time).count()) << ": " << uint2string(_p_min) << "\r";
			if (quit || (std::chrono::duration<double>(now - _record_time).count() > 5 * 60))
			{
				_record_time = now;
				write(quit);
			}
		}
		return !quit;
	}

	bool check_root(const uint64_t b, const uint64_t b_min, const uint64_t p, const int n)
	{
		const size_t i = (b - b_min) / 2;
		if (!_bsieve[i])
		{
			const Mod mod(p);
			if (mod.pow(b, 1 << n) == p - 1)
			{
				_bsieve[i] = true;
				// std::ofstream resFile("res.txt", std::ios::app); if (resFile.is_open()) { resFile << p << " " << b << std::endl; resFile.close(); }
			}
			// May fail if p is not prime
			else if (mod.isprime())
			{
				std::ostringstream ss; ss << "Calculation error (check): p = " << p << ", b = " << b << ".";
				throw std::runtime_error(ss.str());
			}
			else return false;
		}
		return true;
	}

	bool check_roots(const uint64_t b, const uint64_t b_min, const uint64_t b_max, const uint64_t p, const int n)
	{
		for (uint64_t s = b; s <= b_max; s += 2 * p)
		{
			if (s >= b_min)
			{
				if (!check_root(s, b_min, p, n)) return false;
			}
		}
		return true;
	}

	void gen_k()
	{
		const int n = _n;
		const uint64_t k_min = _p_min >> (n + 1);

		static const size_t wheel_prime_count = 8 * 1024;		// largest prime is 84047
		static const uint32_t wheel_sieve_size = 128 * 1024;
		uint32_t wheel_prm[wheel_prime_count];		// 32K
		uint32_t wheel_j_ptr[wheel_prime_count];	// 32K
		bool wheel_sieve[wheel_sieve_size];			// 128K

		wheel_prm[0] = 3; wheel_prm[1] = 5; wheel_prm[2] = 7;
		size_t i = 3;
		for (uint32_t p = 11; i < wheel_prime_count; p += 2)
		{
			const uint32_t s = uint32_t(std::sqrt(double(p))) + 1;
			uint32_t d; for (d = 3; d <= s; d += 2) if (p % d == 0) break;
			if (d > s) { wheel_prm[i] = p; ++i; }
		}

		// std::cout << _wheel_prm[wheel_prime_count - 1] << std::endl;

		uint32_t wheel_j = uint32_t(k_min % wheel_sieve_size);
		uint64_t wheel_k = k_min - wheel_j;

		for (size_t i = 0; i < wheel_prime_count; ++i)
		{
			uint32_t p = wheel_prm[i], j = 0; while ((((wheel_k + j) << (n + 1)) + 1) % p != 0) ++j;
			if (j >= wheel_sieve_size) throw std::runtime_error("Wheel::init failed.");
			wheel_j_ptr[i] = j;
		}

		kArray karray;
		size_t j = 0;

		do
		{
			for (size_t j = 0; j < wheel_sieve_size; ++j) wheel_sieve[j] = false;

			for (size_t i = 0; i < wheel_prime_count; ++i)
			{
				uint32_t p = wheel_prm[i], j = wheel_j_ptr[i]; for (; j < wheel_sieve_size; j += p) wheel_sieve[j] = true;
				wheel_j_ptr[i] = j - wheel_sieve_size;
			}

			do
			{
				if (!wheel_sieve[wheel_j])
				{
					karray.k[j] = wheel_k + wheel_j;
					j = (j + 1) % karray_size;
					if (j == 0) _kqueue.push(karray);
				}
				++wheel_j;
			}
			while (wheel_j < wheel_sieve_size);

			wheel_j = 0;
			wheel_k += wheel_sieve_size;
		}
		while (wheel_k >= wheel_sieve_size);

		if (j > 0)
		{
			for (; j < karray_size; ++j) karray.k[j] = karray.k[0];
			_kqueue.push(karray);
		}
		_kqueue.end();
	}

	void gen_p()
	{
		const int n = _n;

		pArray parray;
		size_t i = 0;

		kArray karray;
		while (_kqueue.pop(karray))
		{
			for (size_t j = 0; j < karray_size; ++j)
			{
				const uint64_t k = karray.k[j];
				if (k == 0) break;

				const uint64_t p = (k << (n + 1)) + 1;
				const MpArith mp(p);

				if (mp.prp())
				{
					uint64_t a = 3, ma = mp.three();
					while (jacobi(a, p) != -1) { ++a; ma = mp.add(ma, mp.one()); if (a == 256) break; }
					if (a < 256)
					{
						const uint64_t c = mp.pow(ma, k), b2 = mp.mul(c, c);
						uint64_t b = mp.toInt(c);

						Mpbb2 & mpbb2 = parray.mpbb2[i];
						mpbb2.p = p; mpbb2.q = mp.q(); mpbb2.b = b; mpbb2.b2 = b2;
						i = (i + 1) % parray_size;
						if (i == 0) _pqueue.push(parray);
					}
				}
			}
		}

		if (i > 0)
		{
			for (; i < parray_size; ++i) parray.mpbb2[i] = parray.mpbb2[0];
			_pqueue.push(parray);
		}
		_pqueue.end();
	}

public:
	void check()
	{
		const int n = _n;
		const uint64_t b_min = _b_min, b_max = _b_max;

		std::cout << "n = " << n << ", b in [" << b_min << ", " << b_max << "]." << std::endl;

		if (!read()) _p_min = (1ull << (n + 1)) + 1;

		if (!init()) return;

		std::thread t_gen_k([=] { gen_k(); }); t_gen_k.detach();
		std::this_thread::sleep_for(std::chrono::milliseconds(200));

		std::thread t_gen_p([=] { gen_p(); }); t_gen_p.detach();
		std::this_thread::sleep_for(std::chrono::milliseconds(200));

		pArray parray;
		while (_pqueue.pop(parray))
		{
			for (size_t j = 0; j < parray_size; j += 4)
			{
				const Mpbb2 & mpbb2_0 = parray.mpbb2[j + 0];
				const Mpbb2 & mpbb2_1 = parray.mpbb2[j + 1];
				const Mpbb2 & mpbb2_2 = parray.mpbb2[j + 2];
				const Mpbb2 & mpbb2_3 = parray.mpbb2[j + 3];
				vu64 vp; vp[0] = mpbb2_0.p; vp[1] = mpbb2_1.p; vp[2] = mpbb2_2.p; vp[3] = mpbb2_3.p;
				vu64 vq; vq[0] = mpbb2_0.q; vq[1] = mpbb2_1.q; vq[2] = mpbb2_2.q; vq[3] = mpbb2_3.q;
				vu64 vb; vb[0] = mpbb2_0.b; vb[1] = mpbb2_1.b; vb[2] = mpbb2_2.b; vb[3] = mpbb2_3.b;
				vu64 vb2; vb2[0] = mpbb2_0.b2; vb2[1] = mpbb2_1.b2; vb2[2] = mpbb2_2.b2; vb2[3] = mpbb2_3.b2;

				MpArith4 mp(vp, vq);

				if (vp[0] <= b_max)
				{
					for (uint64_t i = 0; i < (uint64_t(1) << (n - 1)); ++i)
					{
						for (size_t i = 0; i < 4; ++i)
						{
							const uint64_t beven = (vb[i] % 2 == 0) ? vb[i] : vp[i] - vb[i];	// 0 <= b_even < p
							check_roots(beven, b_min, b_max, vp[i], n);
							check_roots(2 * vp[i] - beven, b_min, b_max, vp[i], n);				// p < 2p - b_even <= 2p
						}

						mp.mul(vb, vb2);
					}
				}
				else
				{
					for (uint64_t i = 0; i < (uint64_t(1) << (n - 1)); ++i)
					{
						vu64 beven; for (size_t i = 0; i < 4; ++i) beven[i] = (vb[i] % 2 == 0) ? vb[i] : vp[i] - vb[i];

						mp.mul(vb, vb2);

						for (size_t i = 0; i < 4; ++i)
						{
							if ((beven[i] >= b_min) && (beven[i] <= b_max)) check_root(beven[i], b_min, vp[i], n);
						}

					}
				}
			}

			_p_min = parray.mpbb2[parray_size - 1].p;
			if (!monitor()) return;
		}

		write(true);
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
