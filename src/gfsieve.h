/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include "ocl.h"
#include "timer.h"
#include "engine.h"
#include "uint96.h"

#include <thread>
#include <mutex>
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
	bool _64bit = false;
	bool _display = false;
	size_t _factorSize = 0;
	uint32_t _n = 0;
	size_t _factorsLoop = 0;
	size_t _savedCount = 0;
	int _log2GlobalWorkSize = 18;
	size_t _localWorkSize = 0;
	timer::time _startTime;
	std::string _extension;
	std::mutex _factor_mutex;
	std::vector<cl_ulong2> _factor;

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

		hFile << "#pragma once" << std::endl << std::endl;
		hFile << "#include <cstdint>" << std::endl << std::endl;

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
	void initEngine(engine & engine, const int log2Global) const
	{
		const size_t globalWorkSize = size_t(1) << log2Global;

		std::stringstream src;
		src << "#define\tlog2GlobalWorkSize\t" << log2Global << std::endl;
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
	void readFactors(engine & engine)
	{
		const size_t factorCount = engine.readFactorCount();
		const double elapsedTime = timer::diffTime(timer::currentTime(), _startTime);

		if (factorCount >= _factorSize) throw std::runtime_error("factor count is too large");

		const std::lock_guard<std::mutex> lock(_factor_mutex);

		const std::string runtime = timer::formatTime(elapsedTime);
		std::cout << factorCount << " factors, time = " << runtime << std::endl;

		if (_savedCount != factorCount)
		{
			_factor.resize(factorCount);
			engine.readFactors(_factor.data(), factorCount);
		}
	}

private:
	void printFactors(const uint64_t cnt)
	{
		const std::lock_guard<std::mutex> lock(_factor_mutex);

		const size_t factorCount = _factor.size();
		if (_savedCount != factorCount)
		{
			const std::string resFilename = std::string("gf") + _extension;
			std::ofstream resFile(resFilename, std::ios::app);
			if (resFile.is_open())
			{
				const uint32_t n = _n, N = uint32_t(1) << n;
				for (size_t i = _savedCount; i < factorCount; ++i)
				{
					const cl_ulong2 & f = _factor[i];
					const uint32_t b = uint32_t(f.s[1] >> 32);
					const uint96 p = uint96(uint32_t(f.s[0]), uint32_t(f.s[0] >> 32), uint32_t(f.s[1]));
					char str[32]; p.get_str(str);

					uint96 r = uint96(b);
					for (size_t j = 0; j < n; ++j) r.square_mod(p);
					r |= 1;
					if (r == p)
					{
						std::ostringstream ss; ss << str << " | " << b << "^" << N << "+1" << std::endl;
						resFile << ss.str();
						if (_display) std::cout << ss.str();
					}
					else
					{
						const uint96 r2 = uint96::pow_mod(2, p, p);
						if (r2 != uint96(2))
						{
							std::ostringstream ss; ss << str << " is not 2-prp";
							throw std::runtime_error(ss.str());
						}
						bool isPrime = true;
						for (uint32_t a = 3; a < 1000; a += 2)
						{
							const uint96 ra = uint96::pow_mod(a, p, p);
							if (ra != uint96(a))
							{
								isPrime = false;
								break;
							}
						}
						if (isPrime)
						{
							std::ostringstream ss; ss << str << " doesn't divide " << b << "^" << N << "+1";
							throw std::runtime_error(ss.str());
						}
					}
				}
				resFile.close();
				_savedCount = factorCount;

				const std::string ctxFilename = std::string("ctx") + _extension;
				std::ofstream ctxFile(ctxFilename);
				if (ctxFile.is_open())
				{
					ctxFile << cnt << " " << _log2GlobalWorkSize << " " << _localWorkSize << std::endl;
					ctxFile.close();
				}
			}
		}
	}

private:
	void saveFactors(engine & engine, const uint64_t cnt, const bool wait)
	{
		readFactors(engine);

		if (wait)
		{
			printFactors(cnt);
		}
		else
		{
			std::thread t( [=] { printFactors(cnt); } );
			t.detach();
		}
	}

private:
	double autoTuning(engine & engine, const uint32_t p_min)
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
				if (_quit) return 0;

				engine.checkPrimes(global, i);
				engine.initFactors(global);
				engine.checkFactors(global, N_2_factors_loop, (local == 8) ? 0 : local);
				engine.clearPrimes();
				engine.readFactorCount();

				const double time = engine.getProfileTime() / double(global);
				if (time < bestTime)
				{
					bestTime = time;
					_log2GlobalWorkSize = log2Global;
					_localWorkSize = (local == 8) ? 0 : local;
				}

				engine.resetProfiles();
			}

			clearEngine(engine);
		}

		const double pTime = bestTime / 1e9 * 1e15 / pow(2.0, double(_n + 1));
		std::ostringstream ss; ss << std::setprecision(3) << 86400 / pTime << " P/day (globalWorkSize = "
			<< (size_t(1) << _log2GlobalWorkSize) << ", localWorkSize = " << _localWorkSize << ")";
		std::cout << ss.str() << std::endl;
		return pTime;
	}

public:
	bool check(engine & engine, const uint32_t n, const uint32_t p_min, const uint32_t p_max, const bool display)
	{
		_64bit = (p_max + 1 <= 9223);	// 2^63 / 10^15
		_display = display;
		_factorSize = (p_min >= 8) ? (size_t(1) << 24) : (size_t(1) << 26);
		_n = n;
		_factorsLoop = size_t(1) << std::min(_n - 1, 10u);
		_savedCount = 0;
		std::stringstream ss; ss << n << "_" << p_min << "_" << p_max << ".txt";
		_extension = ss.str();

		std::cout << (_64bit ? "64" : "96") << "-bit mode" << std::endl;

		uint64_t cnt = 0;
		const std::string ctxFilename = std::string("ctx") + _extension;
		std::ifstream ctxFile(ctxFilename);
		if (ctxFile.is_open())
		{
			ctxFile >> cnt;
			ctxFile >> _log2GlobalWorkSize;
			ctxFile >> _localWorkSize;
			ctxFile.close();
		}

		if (cnt == 0)
		{
			const double pTime = autoTuning(engine, p_min);
			const std::string estimatedTime = timer::formatTime(pTime * (p_max - p_min));
			std::cout << "Estimated time: " << estimatedTime << std::endl;
		}

		engine.setProfiling(false);
		initEngine(engine, _log2GlobalWorkSize);

		const double f = 1e15 / pow(2.0, double(n + 1 + _log2GlobalWorkSize));
		const uint64_t i_min = uint64_t(floor(p_min * f)), i_max = uint64_t(ceil(p_max * f));

		const size_t N_2_factors_loop = (size_t(1) << (_n - 1)) / _factorsLoop;

		const uint64_t i = i_min + cnt;
		uint96 p_min96 = uint96(uint32_t(i), uint32_t(i >> 32));
		p_min96 <<= _log2GlobalWorkSize + _n + 1;
		p_min96 |= 1;

		uint96 p_max96 = uint96(uint32_t(i_max), uint32_t(i_max >> 32));
		p_max96 <<= _log2GlobalWorkSize;
		p_max96 -= 1;
		p_max96 <<= _n + 1;
		p_max96 |= 1;

		char p_min_str[32], p_max_str[32];
		p_min96.get_str(p_min_str); p_max96.get_str(p_max_str);
		std::cout << ((cnt != 0) ? "Resuming from a checkpoint, t" : "T") << "esting n = " << _n << " from " << p_min_str << " to " << p_max_str << std::endl;

		_startTime = timer::currentTime();
		timer::time displayTime = _startTime, recordTime = _startTime;

		const size_t globalWorkSize = size_t(1) << _log2GlobalWorkSize, localWorkSize = _localWorkSize;

		for (uint64_t i = i_min + cnt; i < i_max; ++i)
		{
			if (_quit) break;

			engine.checkPrimes(globalWorkSize, i);
			// const size_t primeCount = engine.readPrimeCount();
			// std::cout << primeCount << " primes" << std::endl;
			engine.initFactors(globalWorkSize);
			engine.checkFactors(globalWorkSize, N_2_factors_loop, localWorkSize);
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
				saveFactors(engine, cnt, false);
			}
		}

		if (cnt > 0)
		{
			std::cout << " terminating...         \r";
			saveFactors(engine, cnt, true);
		}

		// engine.displayProfiles(1);

		clearEngine(engine);

		return true;
	}
};
