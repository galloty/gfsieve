/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <chrono>
#include <string>
#include <sstream>
#include <iomanip>

struct timer
{
	typedef std::chrono::high_resolution_clock::time_point time;

	static time current_time()
	{
		return std::chrono::high_resolution_clock::now();
	}

	static double diff_time(const time & end, const time & start)
	{
		return std::chrono::duration<double>(end - start).count();
	}

	static std::string format_time(const double time)
	{
		uint64_t seconds = static_cast<uint64_t>(time), minutes = seconds / 60, hours = minutes / 60;
		seconds -= minutes * 60; minutes -= hours * 60;

		std::stringstream ss;
		ss << std::setfill('0') << std::setw(2) << hours << ':' << std::setw(2) << minutes << ':' << std::setw(2) << seconds;
		return ss.str();
	}
};
