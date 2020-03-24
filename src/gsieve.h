/*
Copyright 2020, Yves Gallot

gsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include "ocl.h"
#include "timer.h"
#include "engine.h"

#include <thread>
#include <chrono>

#include "ocl/sieve.h"

class gsieve
{
private:
	struct deleter { void operator()(const gsieve * const p) { delete p; } };

public:
	gsieve() {}
	virtual ~gsieve() {}

	static gsieve & getInstance()
	{
		static std::unique_ptr<gsieve, deleter> pInstance(new gsieve());
		return *pInstance;
	}

public:
	void quit() { _quit = true; }

protected:
	volatile bool _quit = false;

private:
	bool readOpenCL(const char * const clFileName, const char * const headerFileName, const char * const varName, std::stringstream & src) const
	{
		std::ifstream clFile(clFileName);
		if (!clFile.is_open()) return false;
		
		// if .cl file exists then generate header file
		std::ofstream hFile(headerFileName, std::ios::binary);	// binary: don't convert line endings to `CRLF` 
		if (!hFile.is_open()) throw std::runtime_error("cannot write openCL header file");

		hFile << "/*" << std::endl;
		hFile << "Copyright 2020, Yves Gallot" << std::endl << std::endl;
		hFile << "gsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it." << std::endl;
		hFile << "Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful." << std::endl;
		hFile << "*/" << std::endl << std::endl;

		hFile << "static const char * const " << varName << " = \\" << std::endl;

		std::string line;
		while (std::getline(clFile, line))
		{
			hFile << "\"";
			for (char c : line)
			{
				if ((c == '\\') || (c == '\"')) hFile << '\\';
				hFile << c;
			}
			hFile << "\\n\" \\" << std::endl;

			src << line << std::endl;
		}
		hFile << "\"\";" << std::endl;

		hFile.close();
		clFile.close();
		return true;
	}

private:
	void uint128_get_str(const __uint128_t & x, char * const str)
	{
		char dgt[32];
		__uint128_t t = x;
		size_t n = 32;
		while (t != 0)
		{
			--n; dgt[n] = '0' + char(t % 10);
			t /= 10;
		}
		for (size_t i = n; i < 32; ++i) str[i - n] = dgt[i];
		str[32 - n] = '\0';
	}

public:
	bool check(const uint32_t n, const uint32_t p_min, const uint32_t p_max, engine & engine)
	{
		const int log2_prime_size = 17;
		const size_t prime_size = size_t(1) << log2_prime_size;
		const size_t factor_size = size_t(1) << 20;

		std::stringstream src;
		src << "#define\tlog2_prime_size\t" << log2_prime_size << std::endl;
		src << "#define\tgfn_n\t" << n << std::endl;
		src << std::endl;

		if (!readOpenCL("ocl/sieve.cl", "src/ocl/sieve.h", "src_ocl_sieve", src)) src << src_ocl_sieve;

		engine.loadProgram(src.str());
		engine.allocMemory(prime_size, factor_size);
		engine.createKernels();

		engine.clearCounters();

		const uint64_t i_min = uint64_t(p_min) << (50 - (n + 1) - log2_prime_size);
		/*const*/ uint64_t i_max = uint64_t(p_max) << (50 - (n + 1) - log2_prime_size);
		const uint64_t i_size = i_max - i_min + 1;

		// std::cout << "i_min = " << i_min << ", i_max = " << i_max << ", " << i_size << std::endl;
		// i_max = i_min + 10;

		const timer::time startTime = timer::currentTime();

		for (uint64_t i = i_min; i <= i_max; ++i)
		{
			engine.checkPrimes(prime_size, i);
			// const size_t prime_count = engine.readPrimeCount();
			// std::cout << prime_count << " primes" << std::endl;
			engine.initFactors(prime_size);
			engine.checkFactors(prime_size, n);
			// const size_t factor_count = engine.readFactorCount();
			// std::cout << factor_count << " factors" << std::endl;
			engine.clearPrimes();
			if (_quit) break;
			// std::cout << (i - i_min) * 100.0 / i_size << std::endl;
		}

		const size_t factor_count = engine.readFactorCount();
		const double elapsed_time = timer::diffTime(timer::currentTime(), startTime) * i_size / (i_max - i_min + 1);
		const std::string runtime = timer::formatTime(elapsed_time);
		std::ostringstream ss; ss << std::setprecision(3) << 86400 / elapsed_time << " Pi/day";

		std::cout << factor_count << " factors, time = " << runtime << ", " << ss.str() << std::endl;

		if (factor_count >= factor_size) throw std::runtime_error("Internal error detected");

		const uint32_t N = uint32_t(1) << n;

		std::vector<cl_ulong2> factor(factor_count);
		engine.readFactors(factor.data(), factor_count);

		engine.releaseKernels();
		engine.releaseMemory();
		engine.clearProgram();

		std::ofstream resFile("res.txt", std::ios::app);
		if (resFile.is_open())
		{
			for (const cl_ulong2 f : factor)
			{
				const uint32_t b = uint32_t(f.s[1] >> 32), f_h = uint32_t(f.s[1]);
				const __uint128_t p = (__uint128_t(f_h) << 64) | f.s[0];
				char str[32]; uint128_get_str(p, str);
				resFile << str << " | " << b << "^" << N << "+1" << std::endl;
			}
			resFile.close();
		}

		return true;
	}
};
