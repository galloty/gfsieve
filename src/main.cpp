/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include "ocl.h"
#include "engine.h"
#include "gfsieve.h"

#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <iostream>

#if defined(_WIN32)
#include <Windows.h>
#else
#include <signal.h>
#endif

class application
{
private:
	struct deleter { void operator()(const application * const p) { delete p; } };

private:
	static void quit(int)
	{
		gfsieve::getInstance().quit();
	}

private:
#if defined(_WIN32)
	static BOOL WINAPI HandlerRoutine(DWORD)
	{
		quit(1);
		return TRUE;
	}
#endif

public:
	application()
	{
#if defined(_WIN32)
		SetConsoleCtrlHandler(HandlerRoutine, TRUE);
#else
		signal(SIGTERM, quit);
		signal(SIGINT, quit);
#endif
	}

	virtual ~application() {}

	static application & getInstance()
	{
		static std::unique_ptr<application, deleter> pInstance(new application());
		return *pInstance;
	}

private:
	static std::string header(const bool nl = false)
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
		ss << "gfsieve 24.9.0 " << sysver << ssc.str() << std::endl;
		ss << "Copyright (c) 2020, Yves Gallot" << std::endl;
		ss << "gfsieve is free source code, under the MIT license." << std::endl;
		if (nl) ss << std::endl;
		return ss.str();
	}

private:
	static std::string usage()
	{
		std::ostringstream ss;
		ss << "Usage: gfsieve <n> <p_min> <p_max> [options]" << std::endl;
		ss << "  <n>                     GFN exponent: b^{2^n} + 1, 15 <= n <= 24" << std::endl;
		ss << "  <p_min>                 start of the p range (unit is P = 10^15)" << std::endl;
		ss << "  <p_max>                 end of the p range (unit is P = 10^15)" << std::endl;
		ss << "  -d <n> or --device <n>  set device number=<n> (default 0)" << std::endl;
		ss << "  -p                      display factors on the screen (default false)" << std::endl;
		ss << "  -v or -V                print the startup banner and immediately exit" << std::endl;
		ss << std::endl;
		return ss.str();
	}

public:
	void run(const std::vector<std::string> & args)
	{
		// if -v or -V then print header to stderr and exit
		for (const std::string & arg : args)
		{
			if ((arg[0] == '-') && ((arg[1] == 'v') || (arg[1] == 'V')))
			{
				std::cerr << header();
				return;
			}
		}

		std::cout << header(true);

		if (args.size() < 3) std::cout << usage();	// print usage, display devices and exit

		platform pfm;
		pfm.displayDevices();

		if (args.size() < 3) return;

		// parse args
		const uint32_t n = (args.size() > 0) ? std::atoi(args[0].c_str()) : 15;	// 23
		const uint32_t p_min = (args.size() > 1) ? std::atoi(args[1].c_str()) : 100000000;	// 9220000, 100000000
		const uint32_t p_max = (args.size() > 2) ? std::atoi(args[2].c_str()) : 100000010;	// 9220010, 100000010
		size_t d = 0;
		bool display = false;
		for (size_t i = 3, size = args.size(); i < size; ++i)
		{
			const std::string & arg = args[i];

			if (arg.substr(0, 2) == "-d")
			{
				const std::string dev = ((arg == "-d") && (i + 1 < size)) ? args[++i] : arg.substr(2);
				d = std::atoi(dev.c_str());
				if (d >= pfm.getDeviceCount()) throw std::runtime_error("invalid device number");
			}

			if (arg == "-p") display = true;
		}

		if (n < 15) throw std::runtime_error("n must be greater than 14");
		if (n > 24) throw std::runtime_error("n must be lesser than 25");
		if (p_min < 1) throw std::runtime_error("p_min must be greater than 0");
		if (p_max <= p_min) throw std::runtime_error("p_max must be greater than p_min");

		gfsieve & sieve = gfsieve::getInstance();

		engine eng(pfm, d);
		sieve.check(eng, n, p_min, p_max, display);
	}
};

int main(int argc, char * argv[])
{
	try
	{
		application & app = application::getInstance();

		std::vector<std::string> args;
		for (int i = 1; i < argc; ++i) args.push_back(argv[i]);
		app.run(args);
	}
	catch (const std::runtime_error & e)
	{
		std::ostringstream ss; ss << std::endl << "error: " << e.what() << "." << std::endl;
		std::cerr << ss.str();
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
