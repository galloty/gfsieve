/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include "ocl.h"
#include "timer.h"
#include "engine.h"

#include <thread>
#include <cmath>

#include <gmp.h>

inline void mpz_set_u64(mpz_t & rop, const uint64_t op)
{
	// _LONG_LONG_LIMB is required
	if (op == 0) mpz_set_ui(rop, 0);
	else { mpz_set_ui(rop, 1); rop->_mp_d[0] = op; }
}

#include "ocl/sieve64.h"
#include "ocl/sieve79.h"

class gfsieve
{
private:
	struct deleter { void operator()(const gfsieve * const p) { delete p; } };

public:
	gfsieve() {}
	virtual ~gfsieve() {}

	static gfsieve & get_instance()
	{
		static std::unique_ptr<gfsieve, deleter> instance(new gfsieve());
		return *instance;
	}

	void quit() { _quit = true; }

protected:
	volatile bool _quit = false;
	bool _64bit = false;
	bool _display = false;
	int _n = 0;
	uint_8 _wheel[8];
	static constexpr size_t _max_factor_size = size_t(1) << 25;	// 512MB
	static constexpr size_t _factors_block = size_t(1) << 14u;
	static constexpr int _log2_block_size = 22;	// => > 285000 primes
	timer::time _start_time;
	std::string _extension;
	std::vector<uint_64_2> _factor;
	std::vector<uint_64> _error;
	static constexpr double _unit = 1e15;

private:
	static bool read_OpenCL(const char * const cl_filename, const char * const header_filename, const char * const varname, std::stringstream & src)
	{
		std::ifstream cl_file(cl_filename);
		if (!cl_file.is_open()) return false;

		// if .cl file exists then generate header file
		std::ofstream h_file(header_filename, std::ios::binary);	// binary: don't convert line endings to `CRLF` 
		if (!h_file.is_open()) throw std::runtime_error("cannot write OpenCL header file");

		h_file << "/*" << std::endl;
		h_file << "Copyright 2020, Yves Gallot" << std::endl << std::endl;
		h_file << "gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it." << std::endl;
		h_file << "Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful." << std::endl;
		h_file << "*/" << std::endl << std::endl;

		h_file << "#pragma once" << std::endl << std::endl;
		h_file << "#include <cstdint>" << std::endl << std::endl;

		h_file << "static const char * const " << varname << " = \\" << std::endl;

		std::string line;
		while (std::getline(cl_file, line))
		{
			h_file << "\"";
			for (char c : line)
			{
				if ((c == '\\') || (c == '\"')) h_file << '\\';
				h_file << c;
			}
			h_file << "\\n\" \\" << std::endl;

			src << line << std::endl;
		}
		h_file << "\"\";" << std::endl;

		h_file.close();
		cl_file.close();
		return true;
	}

	void init_engine(engine & eng, const bool is64) const
	{
		std::stringstream src;
		src << "#define G_N\t" << _n + 1 << std::endl;
		src << "#define LN_BLKSZ\t" << _log2_block_size << std::endl;
		src << "#define FBLK\t" << _factors_block << std::endl;
		src << std::endl;
		src << "typedef	uchar\tuint_8;" << std::endl;
		src << "__constant uint_8 wheel[8] = { ";
		for (size_t i = 0; i < 8; ++i) { src << int(_wheel[i]); if (i != 7) src << ", "; }
		src << " };" << std::endl << std::endl;

		if (_64bit)
		{
			if (!read_OpenCL("ocl/sieve64.cl", "src/ocl/sieve64.h", "src_ocl_sieve64", src)) src << src_ocl_sieve64;
		}
		else
		{
			if (!read_OpenCL("ocl/sieve79.cl", "src/ocl/sieve79.h", "src_ocl_sieve79", src)) src << src_ocl_sieve79;
		}

		eng.loadProgram(src.str());
		// prime_size <= 15/8 * 2 / log(1e15) * block_size = 0.1086 * block_size < 0.25 * block_size
		const size_t prime_size = (size_t(1) << _log2_block_size) / 4;
		eng.alloc_memory(prime_size, _max_factor_size, is64);
		eng.create_kernels();

		// Kronecker symbols (i/j) for odd j <= 255.
		int_8 kro[128 * 256];
		mpz_t zj; mpz_init(zj);
		for (uint32_t j = 3; j < 256; j += 2)
		{
			mpz_set_ui(zj, j);
			for (uint32_t i = 0; i < j; ++i)
			{
				kro[256 * ((j - 3) / 2) + i] = int_8(mpz_ui_kronecker(i, zj));
			}
		}
		mpz_clear(zj);
		eng.write_Kronecker(kro);

		eng.init();
	}

	void clear_engine(engine & eng) const
	{
		eng.release_kernels();
		eng.release_memory();
		eng.clearProgram();
	}

	static bool is_prp(const mpz_t & zp)
	{
		mpz_t zt; mpz_init(zt);
		for (uint32_t a = 3; a < 1000; a += 2)
		{
			mpz_set_ui(zt, a); mpz_powm(zt, zt, zp, zp);
			if (mpz_cmp_ui(zt, a) != 0) { mpz_clear(zt); return false; }
		}
		mpz_clear(zt); return true;
	}

	size_t read_factors(engine & eng, const double dp)
	{
		const size_t factor_count = eng.read_factor_count(), error_count = eng.read_error_count();
		const double elapsed_time = timer::diff_time(timer::current_time(), _start_time);

		if (error_count > 0)
		{
			_error.resize(error_count);
			eng.read_errors(_error.data(), error_count);
			std::cout << error_count << " errors" << std::endl;

			const int n = _n;
			mpz_t zp; mpz_init(zp);
			for (const cl_ulong k : _error)
			{
				mpz_set_u64(zp, k); mpz_mul_2exp(zp, zp, mp_bitcnt_t(n + 1)); mpz_add_ui(zp, zp, 1);
				// if (is_prp(zp))
				{
					std::ostringstream ss; ss << k << " * " << "2^" << n + 1 << " + 1: validation failed";
					throw std::runtime_error(ss.str());
				}
			}
			mpz_clear(zp);
		}

		if (factor_count >= _max_factor_size) throw std::runtime_error("factor count is too large");

		std::cout << factor_count << " factor(s)";
		if (elapsed_time > 10)
		{
			std::cout << ", " << std::max(int(dp * 86400 / elapsed_time), 1) << "P/day, time = " << timer::format_time(elapsed_time);
		}
		else std::cout << "                      ";
		std::cout << std::endl;

		if (factor_count > 0)
		{
			_factor.resize(factor_count);
			eng.read_factors(_factor.data(), factor_count);
		}

		return factor_count;
	}

	bool record_factors()
	{
		bool success = false;
		mpz_t zp, zr; mpz_inits(zp, zr, nullptr);

		const std::string res_filename = std::string("gf") + _extension;
		std::ofstream res_file(res_filename, std::ios::app);
		if (res_file.is_open())
		{
			const int n = _n;
			for (const cl_ulong2 f : _factor)
			{
				const uint32_t b = uint32_t(f.s[1]);
				mpz_set_u64(zp, f.s[0]); mpz_mul_2exp(zp, zp, mp_bitcnt_t(n + 1)); mpz_add_ui(zp, zp, 1);
				char str[32]; mpz_get_str(str, 10, zp);

				mpz_set_ui(zr, b); mpz_powm_ui(zr, zr, uint32_t(1) << n, zp); mpz_add_ui(zr, zr, 1);
				if (mpz_cmp(zr, zp) == 0)
				{
					std::ostringstream ss; ss << str << " | " << b << "^{2^" << n << "}+1" << std::endl;
					res_file << ss.str();
					if (_display) std::cout << ss.str();
				}
				else
				{
					mpz_set_ui(zr, 2); mpz_powm(zr, zr, zp, zp);
					if (mpz_cmp_ui(zr, 2) != 0)
					{
						std::ostringstream ss; ss << str << " is not 2-prp";
						throw std::runtime_error(ss.str());
					}
					if (is_prp(zp))
					{
						std::ostringstream ss; ss << str << " doesn't divide " << b << "^{2^" << n << "}+1";
						throw std::runtime_error(ss.str());
					}
				}
			}
			res_file.close();
			success = true;
		}

		mpz_clears(zp, zr, nullptr);
		return success;
	}

	void save_factors(engine & eng, const uint64_t cnt, const double dp)
	{
		if (read_factors(eng, dp) > 0)
		{
			if (record_factors())
			{
				const std::string ctx_filename = std::string("ctx") + _extension;
				std::ofstream ctx_file(ctx_filename);
				if (ctx_file.is_open())
				{
					ctx_file << cnt << std::endl;
					ctx_file.close();
				}
				eng.clear_factor_count();
			}
		}
	}

	uint_64 get_k(const uint_64 i) const
	{
		const uint_64 j = i << _log2_block_size;
		return 15 * (j / 8) + _wheel[j % 8];
	}

public:
	bool check(engine & eng, const int n, const uint32_t p_min, const uint32_t p_max, const bool display)
	{
		_64bit = (p_max <= size_t(-1) / _unit);
		_display = display;
		_n = n;
		std::stringstream ss; ss << n << "_" << p_min << "_" << p_max << ".txt";
		_extension = ss.str();

		// gcd(p mod 15, 15) = 1: 8 solutions
		size_t i = 0;
		for (uint_64 k = 0; k < 15; ++k)
		{
			const uint_64 p = (k << (n + 1)) | 1;
			if ((p % 3 == 0) || (p % 5 == 0)) continue;
			_wheel[i] = uint_8(k);
			++i;
		}
		if (i != 8) throw std::runtime_error("wheel count != 8");
		std::cout << (_64bit ? 64 : 79) << "-bit mode" << std::endl;

		uint64_t cnt = 0;
		const std::string ctxFilename = std::string("ctx") + _extension;
		std::ifstream ctx_file(ctxFilename);
		if (ctx_file.is_open())
		{
			ctx_file >> cnt;
			ctx_file.close();
		}

		eng.setProfiling(false);
		init_engine(eng, _64bit);

		const size_t N_2_factors_block = (size_t(1) << (n - 1)) / _factors_block;

		const double f = _unit * 8.0 / 15 / std::pow(2.0, double(_log2_block_size + n + 1));
		const uint_64 i_min = uint_64(floor(p_min * f)), i_max = uint_64(ceil(p_max * f));
		const uint_64 i_start = i_min + cnt;

		mpz_t zp_min, zp_max; mpz_inits(zp_min, zp_max, nullptr);

		mpz_set_u64(zp_min, get_k(i_start)); mpz_mul_2exp(zp_min, zp_min, mp_bitcnt_t(n + 1)); mpz_add_ui(zp_min, zp_min, 1);
		mpz_set_u64(zp_max, get_k(i_max)); mpz_mul_2exp(zp_max, zp_max, mp_bitcnt_t(n + 1)); mpz_add_ui(zp_max, zp_max, 1);

		char p_min_str[32], p_max_str[32]; mpz_get_str(p_min_str, 10, zp_min); mpz_get_str(p_max_str, 10, zp_max);
		std::cout << ((cnt != 0) ? "Resuming from a checkpoint, t" : "T") << "esting n = " << n << " from " << p_min_str << " to " << p_max_str << std::endl;

		mpz_clears(zp_min, zp_max, nullptr);

		// std::cout << "For i = " << i_start << " to " << i_max - 1 << std::endl;
		const double log_p_min = std::log(get_k(i_start)) + (n + 1) * std::log(2), log_p_max = std::log(get_k(i_max)) + (n + 1) * std::log(2);
		const double count = (std::log(log_p_max) - std::log(log_p_min)) * 1e9;
		std::cout << "Expected number of factors : ";
		if (count >= 1000) std::cout << uint_64(count); else std::cout << count;
		std::cout << std::endl;
		if (count > 0.9 * _max_factor_size) throw std::runtime_error("range is too large");

		_start_time = timer::current_time();
		timer::time display_time = _start_time, record_time = _start_time;

		const size_t block_size = size_t(1) << _log2_block_size, prime_size = block_size / 4;

		for (uint_64 i = i_start; i < i_max; ++i)
		{
			if (_quit) break;

			eng.generate_primes(block_size, i);
			// const size_t prime_count = eng.read_prime_count();
			// std::cout << i << ": " << prime_count << " primes" << std::endl;
			eng.init_factors(prime_size);
			eng.check_factors(prime_size, N_2_factors_block);
			// const size_t factor_count = eng.read_factor_count();
			// std::cout << factor_count << " factor(s)" << std::endl;
			eng.clear_prime_count();
			++cnt;

			const timer::time current_time = timer::current_time();

			if (timer::diff_time(current_time, display_time) > 1)
			{
				display_time = current_time;
				std::ostringstream ss; ss << std::setprecision(3) << " " << cnt * 100.0 / (i_max - i_min) << "% done    \r";
				std::cout << ss.str();
			}

			if (timer::diff_time(current_time, record_time) > 300)
			{
				record_time = current_time;
				save_factors(eng, cnt, (i - i_start) / f);
			}
		}

		if (cnt > 0)
		{
			std::cout << " terminating...         \r";
			save_factors(eng, cnt, (i_max - i_start) / f);
		}

		// eng.displayProfiles(1);

		clear_engine(eng);

		return true;
	}
};
