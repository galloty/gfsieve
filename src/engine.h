/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include "ocl.h"

typedef cl_uchar	uint_8;
typedef cl_char		int_8;
typedef cl_uint		uint_32;
typedef cl_ulong	uint_64;
typedef cl_ulong2	uint_64_2;

class engine : public device
{
private:
	cl_mem _kro_vector = nullptr;
	cl_mem _k_vector = nullptr, _q_vector = nullptr, _ext_vector = nullptr, _c_vector = nullptr, _cn_vector = nullptr;
	cl_mem _factor_vector = nullptr, _error_vector = nullptr;
	cl_mem _prime_count = nullptr, _factor_count = nullptr, _error_count = nullptr;
	cl_kernel _generate_primes = nullptr, _init_factors = nullptr, _check_factors = nullptr, _clear = nullptr;

public:
	engine(const platform & pfm, const size_t d) : device(pfm, d, true) {}
	virtual ~engine() {}

	void alloc_memory(const size_t prime_size, const size_t factor_size, const bool is64)
	{
#if defined(ocl_debug)
		std::cerr << "Alloc gpu memory." << std::endl;
#endif
		const size_t type_size = is64 ? sizeof(uint_64) : sizeof(uint_64_2);
		_kro_vector = _createBuffer(CL_MEM_READ_ONLY, 128 * 256 * sizeof(int_8));
		_k_vector = _createBuffer(CL_MEM_READ_WRITE, sizeof(uint_64) * prime_size);
		_q_vector = _createBuffer(CL_MEM_READ_WRITE, type_size * prime_size);
		_ext_vector = _createBuffer(CL_MEM_READ_WRITE, type_size * prime_size);
		_c_vector = _createBuffer(CL_MEM_READ_WRITE, type_size * prime_size);
		_cn_vector = _createBuffer(CL_MEM_READ_WRITE, type_size * prime_size);
		_factor_vector = _createBuffer(CL_MEM_READ_WRITE, sizeof(uint_64_2) * factor_size);
		_error_vector = _createBuffer(CL_MEM_READ_WRITE, sizeof(uint_64) * 1048576);
		_prime_count = _createBuffer(CL_MEM_READ_WRITE, sizeof(uint_32));
		_factor_count = _createBuffer(CL_MEM_READ_WRITE, sizeof(uint_32));
		_error_count = _createBuffer(CL_MEM_READ_WRITE, sizeof(uint_32));
	}

	void release_memory()
	{
#if defined(ocl_debug)
		std::cerr << "Free gpu memory." << std::endl;
#endif
		_releaseBuffer(_kro_vector);
		_releaseBuffer(_k_vector);
		_releaseBuffer(_q_vector);
		_releaseBuffer(_ext_vector);
		_releaseBuffer(_c_vector);
		_releaseBuffer(_cn_vector);
		_releaseBuffer(_factor_vector);
		_releaseBuffer(_error_vector);
		_releaseBuffer(_prime_count);
		_releaseBuffer(_factor_count);
		_releaseBuffer(_error_count);
	}

	void create_kernels()
	{
#if defined(ocl_debug)
		std::cerr << "Create ocl kernels." << std::endl;
#endif
		_generate_primes = _createKernel("generate_primes");
		_setKernelArg(_generate_primes, 0, sizeof(cl_mem), &_prime_count);
		_setKernelArg(_generate_primes, 1, sizeof(cl_mem), &_k_vector);
		_setKernelArg(_generate_primes, 2, sizeof(cl_mem), &_q_vector);
		_setKernelArg(_generate_primes, 3, sizeof(cl_mem), &_ext_vector);

		_init_factors = _createKernel("init_factors");
		_setKernelArg(_init_factors, 0, sizeof(cl_mem), &_prime_count);
		_setKernelArg(_init_factors, 1, sizeof(cl_mem), &_k_vector);
		_setKernelArg(_init_factors, 2, sizeof(cl_mem), &_q_vector);
		_setKernelArg(_init_factors, 3, sizeof(cl_mem), &_ext_vector);
		_setKernelArg(_init_factors, 4, sizeof(cl_mem), &_kro_vector);
		_setKernelArg(_init_factors, 5, sizeof(cl_mem), &_c_vector);
		_setKernelArg(_init_factors, 6, sizeof(cl_mem), &_cn_vector);

		_check_factors = _createKernel("check_factors");
		_setKernelArg(_check_factors, 0, sizeof(cl_mem), &_prime_count);
		_setKernelArg(_check_factors, 1, sizeof(cl_mem), &_k_vector);
		_setKernelArg(_check_factors, 2, sizeof(cl_mem), &_q_vector);
		_setKernelArg(_check_factors, 3, sizeof(cl_mem), &_c_vector);
		_setKernelArg(_check_factors, 4, sizeof(cl_mem), &_ext_vector);
		_setKernelArg(_check_factors, 5, sizeof(cl_mem), &_cn_vector);
		_setKernelArg(_check_factors, 6, sizeof(cl_mem), &_factor_count);
		_setKernelArg(_check_factors, 7, sizeof(cl_mem), &_factor_vector);
		_setKernelArg(_check_factors, 8, sizeof(cl_mem), &_error_count);
		_setKernelArg(_check_factors, 9, sizeof(cl_mem), &_error_vector);

		_clear = _createKernel("clear");
	}

	void release_kernels()
	{
#if defined(ocl_debug)
		std::cerr << "Release ocl kernels." << std::endl;
#endif
		_releaseKernel(_generate_primes);
		_releaseKernel(_init_factors);
		_releaseKernel(_check_factors);
		_releaseKernel(_clear);
	}

	void init()
	{
		const uint_32 zero = 0;
		_writeBuffer(_prime_count, &zero, sizeof(uint_32));
		_writeBuffer(_factor_count, &zero, sizeof(uint_32));
		_writeBuffer(_error_count, &zero, sizeof(uint_32));
	}

	uint_32 read_prime_count() { uint_32 count; _readBuffer(_prime_count, &count, sizeof(uint_32)); return count; }
	uint_32 read_factor_count() { uint_32 count; _readBuffer(_factor_count, &count, sizeof(uint_32)); return count; }
	uint_32 read_error_count() { uint_32 count; _readBuffer(_error_count, &count, sizeof(uint_32)); return count; }

	void read_factors(uint_64_2 * const ptr, const size_t count) { if (count > 0) _readBuffer(_factor_vector, ptr, sizeof(uint_64_2) * count); }
	void read_errors(uint_64 * const ptr, const size_t count) { if (count > 0) _readBuffer(_error_vector, ptr, sizeof(uint_64) * count); }

	void write_Kronecker(const int_8 * const data) { _writeBuffer(_kro_vector, data, 128 * 256 * sizeof(int_8)); }

	void generate_primes(const size_t size, const uint_64 i)
	{
		_setKernelArg(_generate_primes, 4, sizeof(uint_64), &i);
		_executeKernel(_generate_primes, size);
	}

	void init_factors(const size_t size) { _executeKernel(_init_factors, size); }

	void check_factors(const size_t size, const size_t count)
	{
		cl_char last = cl_char(0);
		_setKernelArg(_check_factors, 10, sizeof(cl_char), &last);
		for (size_t i = 0; i < count - 1; ++i) _executeKernel(_check_factors, size);
		last = cl_char(1);
		_setKernelArg(_check_factors, 10, sizeof(cl_char), &last);
		_executeKernel(_check_factors, size);
	}

	void clear_prime_count()
	{
		_setKernelArg(_clear, 0, sizeof(cl_mem), &_prime_count);
		_executeKernel(_clear, 1);
	}

	void clear_factor_count()
	{
		_setKernelArg(_clear, 0, sizeof(cl_mem), &_factor_count);
		_executeKernel(_clear, 1);
		_setKernelArg(_clear, 0, sizeof(cl_mem), &_error_count);
		_executeKernel(_clear, 1);
	}
};
