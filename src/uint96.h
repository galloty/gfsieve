/*
Copyright 2020, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>

inline uint32_t lo32(const uint64_t n) { return uint32_t(n); }
inline uint32_t hi32(const uint64_t n) { return uint32_t(n >> 32); }
inline uint64_t lo(const uint64_t n) { return uint64_t(lo32(n)); }
inline uint64_t hi(const uint64_t n) { return uint64_t(hi32(n)); }

class uint96
{
private:
	uint32_t _a[3];

	struct uint192
	{
		uint32_t _a[6];

		uint192() {}
		virtual ~uint192() {}
		uint192(const uint32_t n) { _a[0] = n; _a[1] = 0; _a[2] = 0; _a[3] = 0; _a[4] = 0; _a[5] = 0; }
		uint192(const uint96 & n) { _a[0] = n._a[0]; _a[1] = n._a[1]; _a[2] = n._a[2]; _a[3] = 0; _a[4] = 0; _a[5] = 0; }

		operator uint96() const { return uint96(_a[0], _a[1], _a[2]); }

		void dup()
		{
			for (size_t i = 5; i > 0; --i) _a[i] = (_a[i] << 1) | (_a[i - 1] >> 31);
			_a[0] <<= 1;
		}

		void half()
		{
			for (size_t i = 0; i < 5; ++i) _a[i] = (_a[i] >> 1) | (_a[i + 1] << 31);
			_a[5] >>= 1;
		}

		bool operator <(const uint192 & rhs) const
		{
			for (size_t i = 6; i > 0; --i)
			{
				if (_a[i - 1] < rhs._a[i - 1]) return true;
				if (_a[i - 1] > rhs._a[i - 1]) return false;
			}
			return false;
		}

		bool operator >=(const uint192 & rhs) const
		{
			for (size_t i = 6; i > 0; --i)
			{
				if (_a[i - 1] > rhs._a[i - 1]) return true;
				if (_a[i - 1] < rhs._a[i - 1]) return false;
			}
			return true;
		}

		uint192 & operator +=(const uint192 & rhs)
		{
			uint64_t l = 0;
			for (size_t i = 0; i < 6; ++i)
			{
				l += uint64_t(_a[i]) + uint64_t(rhs._a[i]);
				_a[i] = uint32_t(l);
				l >>= 32;
			}
			return *this;
		}

		uint192 & operator -=(const uint192 & rhs)
		{
			int64_t l = 0;
			for (size_t i = 0; i < 6; ++i)
			{
				l += int64_t(_a[i]) - int64_t(rhs._a[i]);
				_a[i] = uint32_t(l);
				l >>= 32;
			}
			return *this;
		}

		void square96()
		{
			const uint64_t a00 = _a[0] * uint64_t(_a[0]), a01 = _a[0] * uint64_t(_a[1]), a02 = _a[0] * uint64_t(_a[2]);
			const uint64_t a11 = _a[1] * uint64_t(_a[1]), a12 = _a[1] * uint64_t(_a[2]), a22 = _a[2] * uint64_t(_a[2]);
			const uint64_t s1 = hi(a00) + 2 * lo(a01);
			const uint64_t s2 = lo(a11) + 2 * (hi(a01) + lo(a02)) + hi(s1);
			const uint64_t s3 = hi(a11) + 2 * (hi(a02) + lo(a12)) + hi(s2);
			const uint64_t s4 = a22 + 2 * hi(a12) + hi(s3);
			_a[0] = lo32(a00); _a[1] = lo32(s1); _a[2] = lo32(s2);
			_a[3] = lo32(s3); _a[4] = lo32(s4); _a[5] = hi32(s4);
		}

		void mul96(const uint32_t n)
		{
			const uint64_t a00 = _a[0] * uint64_t(n), a01 = _a[1] * uint64_t(n), a02 = _a[2] * uint64_t(n);
			const uint64_t s1 = hi(a00) + lo(a01);
			const uint64_t s2 = a02 + hi(a01) + hi(s1);
			_a[0] = lo32(a00); _a[1] = lo32(s1); _a[2] = lo32(s2); _a[3] = hi32(s2);
		}

		void rem(const uint96 & m)
		{
			uint192 num = uint192(m);
			int c = 0;

			while (num < *this)
			{
				num.dup();
				++c;
			}

			while (c != 0)
			{
				if (*this >= num) { *this -= num; }
				num.half();
				--c;
			}

			if (*this >= num) { *this -= num; }
		}
	};

public:
	uint96() {}
	uint96(const uint32_t a0, const uint32_t a1 = 0, const uint32_t a2 = 0) { _a[0] = a0; _a[1] = a1; _a[2] = a2; }
	virtual ~uint96() {}

	bool operator ==(const uint96 & rhs) const
	{
		for (size_t i = 0; i < 3; ++i) if (_a[i] != rhs._a[i]) return false;
		return true;
	}

	bool operator !=(const uint96 & rhs) const { return !(*this == rhs); }

	bool bit(const int s) const
	{
		if (s < 32) return ((_a[0] & (uint32_t(1) << s)) != 0);
		if (s < 64) return ((_a[1] & (uint32_t(1) << (s - 32))) != 0);
		return ((_a[2] & (uint32_t(1) << (s - 64))) != 0);
	}

	uint96 & operator |=(const uint32_t n) { _a[0] |= n; return *this; }

	uint96 & operator +=(const uint32_t n)
	{
		uint64_t l = _a[0] + uint64_t(n);
		_a[0] = uint32_t(l);
		l = _a[1] + (l >> 32);
		_a[1] = uint32_t(l);
		l = _a[2] + (l >> 32);
		_a[2] = uint32_t(l);
		return *this;
	}

	uint96 & operator -=(const uint32_t n)
	{
		int64_t l = -int64_t(n);
		for (size_t i = 0; i < 3; ++i)
		{
			l += int64_t(_a[i]);
			_a[i] = uint32_t(l);
			l >>= 32;
		}
		return *this;
	}

	uint96 & operator <<=(const int s)
	{
		if (s < 32)
		{
			_a[2] = (_a[2] << s) | (_a[1] >> (32 - s));
			_a[1] = (_a[1] << s) | (_a[0] >> (32 - s));
			_a[0] <<= s;
		}
		else if (s == 32)
		{
			_a[2] = _a[1];
			_a[1] = _a[0];
			_a[0] = 0;
		}
		else
		{
			_a[2] = (_a[1] << (s - 32)) | (_a[0] >> (64 - s));
			_a[1] = _a[0] << (s - 32);
			_a[0] = 0;
		}
		return *this;
	}

	void square_mod(const uint96 & m)
	{
		uint192 l = uint192(*this);
		l.square96();
		l.rem(m);
		*this = uint96(l);
	}

	void mul_mod(const uint32_t n, const uint96 & m)
	{
		uint192 l = uint192(*this);
		l.mul96(n);
		l.rem(m);
		*this = uint96(l);
	}

	void get_str(char * const str) const
	{
		char dgt[32];
		uint64_t l48 = _a[0] | (uint64_t(_a[1] & ((uint32_t(1) << 16) - 1)) << 32), h48 = (uint64_t(_a[2]) << 16) | (_a[1] >> 16);
		size_t n = 32;
		while (h48 != 0)
		{
			l48 |= (h48 % 10) << 48;
			h48 /= 10;
			--n; dgt[n] = '0' + char(l48 % 10);
			l48 /= 10;
		}
		while (l48 != 0)
		{
			--n; dgt[n] = '0' + char(l48 % 10);
			l48 /= 10;
		}
		for (size_t i = n; i < 32; ++i) str[i - n] = dgt[i];
		str[32 - n] = '\0';
	}

	static uint96 pow_mod(const uint32_t a, const uint96 & e, const uint96 & m)
	{
		bool s = false;
		uint96 r = uint96(a);
		for (int b = 95; b >= 0; --b)
		{
			if (s) r.square_mod(m);

			if (e.bit(b))
			{
				if (s) r.mul_mod(a, m);
				s = true;
			}
		}
		return r;
	}
};
