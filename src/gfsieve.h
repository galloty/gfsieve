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
#include <chrono>
#include <cmath>

#include "ocl/sieve.h"

class gfsieve
{
private:
	struct deleter { void operator()(const gfsieve * const p) { delete p; } };

public:
	gfsieve() {}
	virtual ~gfsieve() {}

	static gfsieve & getInstance()
	{
		static std::unique_ptr<gfsieve, deleter> pInstance(new gfsieve());
		return *pInstance;
	}

public:
	void quit() { _quit = true; }

protected:
	volatile bool _quit = false;

private:
	static bool readOpenCL(const char * const clFileName, const char * const headerFileName, const char * const varName, std::stringstream & src)
	{
		std::ifstream clFile(clFileName);
		if (!clFile.is_open()) return false;
		
		// if .cl file exists then generate header file
		std::ofstream hFile(headerFileName, std::ios::binary);	// binary: don't convert line endings to `CRLF` 
		if (!hFile.is_open()) throw std::runtime_error("cannot write openCL header file");

		hFile << "/*" << std::endl;
		hFile << "Copyright 2020, Yves Gallot" << std::endl << std::endl;
		hFile << "gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it." << std::endl;
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
	static void uint96_get_str(const uint64_t x_l, const uint32_t x_h, char * const str)
	{
		char dgt[32];
		uint64_t l = x_l & ((uint64_t(1) << 48) - 1), h = (uint64_t(x_h) << 16) | (x_l >> 48);
		size_t n = 32;
		while ((l != 0) || (h != 0))
		{
			l |= (h % 10) << 48;
			h /= 10;
			--n; dgt[n] = '0' + char(l % 10);
			l /= 10;
		}
		for (size_t i = n; i < 32; ++i) str[i - n] = dgt[i];
		str[32 - n] = '\0';
	}

public:
	bool check(const uint32_t n, const uint32_t p_min, const uint32_t p_max, engine & engine) const
	{
		const int log2_prime_size = 18;
		const size_t prime_size = size_t(1) << log2_prime_size;
		const size_t factor_size = size_t(1) << 24;

		std::stringstream src;
		src << "#define\tlog2_prime_size\t" << log2_prime_size << std::endl;
		src << "#define\tgfn_n\t" << n << std::endl;
		src << std::endl;

		if (!readOpenCL("ocl/sieve.cl", "src/ocl/sieve.h", "src_ocl_sieve", src)) src << src_ocl_sieve;

		engine.loadProgram(src.str());
		engine.allocMemory(prime_size, factor_size);
		engine.createKernels();

		engine.clearCounters();

		// engine.setProfiling(true);

		const double f = 1e15 / pow(2.0, double(n + 1 + log2_prime_size));
		const uint64_t i_min = uint64_t(floor(p_min * f)), i_max = uint64_t(ceil(p_max * f));

		const timer::time startTime = timer::currentTime();
		timer::time displayTime = startTime;

		uint64_t cnt = 0;
		for (uint64_t i = i_min; i < i_max; ++i)
		{
			engine.checkPrimes(prime_size, i);
			// const size_t prime_count = engine.readPrimeCount();
			// std::cout << prime_count << " primes" << std::endl;
			engine.initFactors(prime_size);
			engine.checkFactors(prime_size, n);
			// const size_t factor_count = engine.readFactorCount();
			// std::cout << factor_count << " factors" << std::endl;
			engine.clearPrimes();
			++cnt;
			if (_quit) break;
			if (timer::diffTime(timer::currentTime(), displayTime) > 1)
			{
				displayTime = timer::currentTime();
				std::ostringstream ss; ss << std::setprecision(3) << " " << cnt * 100.0 / (i_max - i_min) << "% done    \r";
				std::cout << ss.str();
			}
		}

		std::cout << " Terminating...         \r";
		const size_t factor_count = engine.readFactorCount();
		const double elapsed_time = timer::diffTime(timer::currentTime(), startTime);
		const std::string runtime = timer::formatTime(elapsed_time);
		std::ostringstream ss; ss << std::setprecision(3) << 86400 / (elapsed_time * (i_max - i_min) / cnt) << " P/day";

		std::cout << factor_count << " factors, time = " << runtime << ", " << ss.str() << std::endl;

		if (factor_count >= factor_size) throw std::runtime_error("Internal error detected");

		const uint32_t N = uint32_t(1) << n;

		std::vector<cl_ulong2> factor(factor_count);
		engine.readFactors(factor.data(), factor_count);

		// engine.displayProfiles(cnt);

		engine.releaseKernels();
		engine.releaseMemory();
		engine.clearProgram();

		std::ofstream resFile("res.txt", std::ios::app);
		if (resFile.is_open())
		{
			for (const cl_ulong2 f : factor)
			{
				const uint32_t b = uint32_t(f.s[1] >> 32);
				char str[32]; uint96_get_str(f.s[0], uint32_t(f.s[1]), str);
				resFile << str << " | " << b << "^" << N << "+1" << std::endl;
			}
			resFile.close();
		}

		return true;
	}
};
