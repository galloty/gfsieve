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
#include "ocl/sieve64.h"

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
	const size_t _factorSize = size_t(1) << 24;
	bool _64bit = false;
	uint32_t _n = 0;
	size_t _factorsLoop = 0;
	size_t _savedCount = 0;
	timer::time _startTime;
	std::string _outFilename;

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
	void initEngine(engine & engine, const int log2GlobalWorkSize) const
	{
		const size_t globalWorkSize = size_t(1) << log2GlobalWorkSize;

		std::stringstream src;
		src << "#define\tlog2GlobalWorkSize\t" << log2GlobalWorkSize << std::endl;
		src << "#define\tgfn_n\t" << _n << std::endl;
		src << "#define\tfactors_loop\t" << _factorsLoop << std::endl;
		src << std::endl;

		if (_64bit)
		{
			if (!readOpenCL("ocl/sieve64.cl", "src/ocl/sieve64.h", "src_ocl_sieve64", src)) src << src_ocl_sieve64;
		}
		else
		{
			if (!readOpenCL("ocl/sieve.cl", "src/ocl/sieve.h", "src_ocl_sieve", src)) src << src_ocl_sieve;
		}

		engine.loadProgram(src.str());
		engine.allocMemory(globalWorkSize, _factorSize);
		engine.createKernels();

		engine.clearCounters();
	}

private:
	void clearEngine(engine & engine) const
	{
		engine.releaseKernels();
		engine.releaseMemory();
		engine.clearProgram();
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

private:
	void saveFactors(engine & engine, const double percentDone)
	{
		const size_t factorCount = engine.readFactorCount();
		const double elapsedTime = timer::diffTime(timer::currentTime(), _startTime);

		if (factorCount >= _factorSize) throw std::runtime_error("Internal error detected");

		if (_savedCount != factorCount)
		{
			const uint32_t N = uint32_t(1) << _n;

			std::vector<cl_ulong2> factor(factorCount);
			engine.readFactors(factor.data(), factorCount);

			std::ofstream resFile(_outFilename, std::ios::app);
			if (resFile.is_open())
			{
				for (size_t i = _savedCount; i < factorCount; ++i)
				{
					const cl_ulong2 & f = factor[i];
					const uint32_t b = uint32_t(f.s[1] >> 32);
					char str[32]; uint96_get_str(f.s[0], uint32_t(f.s[1]), str);
					std::ostringstream ss; ss << str << " | " << b << "^" << N << "+1" << std::endl;
					resFile << ss.str();
					std::cout << ss.str();
				}
				resFile.close();
				_savedCount = factorCount;
			}
		}

		std::ostringstream ss; ss << std::setprecision(3) << 86400 / (elapsedTime * percentDone) << " P/day";
		const std::string runtime = timer::formatTime(elapsedTime);
		std::cout << factorCount << " factors, time = " << runtime << ", " << ss.str() << std::endl;
	}

private:
	void autoTuning(engine & engine, const uint32_t p_min, int & log2GlobalWorkSize, size_t & localWorkSize)
	{
		std::cout << " auto-tuning...\r";

		const size_t maxWorkGroupSize = engine.getMaxWorkGroupSize();
		const size_t N_2_factors_loop = (size_t(1) << (_n - 1)) / _factorsLoop;

		double bestTime = 1e100;
		for (int log2Global = 17; log2Global <= 21; ++log2Global)
		{
			const size_t global = size_t(1) << log2Global;

			engine.setProfiling(true);
			initEngine(engine, log2Global);
			
			const double f = 1e15 / pow(2.0, double(_n + 1 + log2Global));
			const uint64_t i = uint64_t(floor(p_min * f));

			for (size_t local = 8; local <= maxWorkGroupSize; local *= 2)
			{
				if (_quit) return;

				engine.checkPrimes(global, i);
				engine.initFactors(global);
				engine.checkFactors(global, N_2_factors_loop, (local == 8) ? 0 : local);
				engine.clearPrimes();
				engine.readFactorCount();

				const double time = engine.getProfileTime() / double(global);
				if (time < bestTime)
				{
					bestTime = time;
					log2GlobalWorkSize = log2Global;
					localWorkSize = (local == 8) ? 0 : local;
				}

				engine.resetProfiles();
			}

			clearEngine(engine);
		}
	}

public:
	bool check(engine & engine, const uint32_t n, const uint32_t p_min, const uint32_t p_max)
	{
		_64bit = (p_max + 1 <= 9223);	// 2^63 / 10^15
		_n = n;
		_factorsLoop = size_t(1) << std::min(_n - 1, 10u);
		_savedCount = 0;
		std::stringstream ss; ss << "gf" << n << "_" << p_min << "_" << p_max << ".txt";
		_outFilename = ss.str();

		std::cout << (_64bit ? "64" : "96") << "-bit mode" << std::endl;

		int log2GlobalWorkSize = 18;
		size_t localWorkSize = 0;
		autoTuning(engine, p_min, log2GlobalWorkSize, localWorkSize);

		const size_t globalWorkSize = size_t(1) << log2GlobalWorkSize;

		std::cout << "globalWorkSize = " << globalWorkSize << ", localWorkSize = " << localWorkSize << std::endl;

		engine.setProfiling(false);
		initEngine(engine, log2GlobalWorkSize);

		const double f = 1e15 / pow(2.0, double(n + 1 + log2GlobalWorkSize));
		const uint64_t i_min = uint64_t(floor(p_min * f)), i_max = uint64_t(ceil(p_max * f));

		_startTime = timer::currentTime();
		timer::time displayTime = _startTime, recordTime = _startTime;

		const size_t N_2_factors_loop = (size_t(1) << (_n - 1)) / _factorsLoop;

		const double total = double(i_max - i_min) / (p_max - p_min);
		uint64_t cnt = 0;
		for (uint64_t i = i_min; i < i_max; ++i)
		{
			if (_quit) break;

			engine.checkPrimes(globalWorkSize, i);
			// const size_t primeCount = engine.readPrimeCount();
			// std::cout << primeCount << " primes" << std::endl;
			engine.initFactors(globalWorkSize);
			engine.checkFactors(globalWorkSize, N_2_factors_loop, 0);
			// const size_t factorCount = engine.readFactorCount();
			// std::cout << factorCount << " factors" << std::endl;
			engine.clearPrimes();
			++cnt;

			const timer::time currentTime = timer::currentTime();

			if (timer::diffTime(currentTime, displayTime) > 1)
			{
				displayTime = currentTime;
				std::ostringstream ss; ss << std::setprecision(3) << " " << cnt * 100.0 / (i_max - i_min) << "% done    \r";
				std::cout << ss.str();
			}

			if (timer::diffTime(currentTime, recordTime) > 300)
			{
				recordTime = currentTime;
				saveFactors(engine, total / cnt);
			}
		}

		if (cnt > 0)
		{
			std::cout << " terminating...         \r";
			saveFactors(engine, total / cnt);
		}

		// engine.displayProfiles(1);

		clearEngine(engine);

		return true;
	}
};
