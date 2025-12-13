/*
Copyright 2025, Yves Gallot

gfsieve is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

\\ 2^{(p - 1)/2} ?= +/-1 mod p
prp(p)=
{
	t = Mod(2, p)^((p - 1)/2);
	return ((t == 1) || (t == -1));
}

main()=
{
	\\ Kronecker symbols (i/j) for odd j <= 255.
	kro_vector = vectorsmall(128 * 256);
	forstep (j = 3, 255, 2,
		for (i = 0, j - 1,
			kro_vector[1 + 256 * ((j - 3) / 2) + i] = kronecker(i, j);
		);
	);

	n = 23; g_n = n + 1;
	
	\\ gcd(p mod 15, 15) = 1: 8 solutions
	wheel = vectorsmall(8);
	i = 0;
	for (k = 0, 15 - 1,
		p = (k << g_n) + 1;
		if ((p % 3 != 0) && (p % 5 != 0),
			wheel[1 + i] = k;
			i++;
		);
	);
	if (i != 8, print("Error: wheel count."); return;);
	
	log2_global_worksize = 22;

	p_min = 10; p_max = 11; unit = 1e15;

	f = unit * 8.0 / 15 / 2^(n + 1 + log2_global_worksize);
	i_min = floor(p_min * f); i_max = ceil(p_max * f);

	j_min = i_min << log2_global_worksize; j_max = i_max << log2_global_worksize;
	k_min = 15 * (j_min \ 8) + wheel[1 + j_min % 8]; k_max = 15 * (j_max \ 8) + wheel[1 + j_max % 8];
	zp_min = (k_min << g_n) + 1; zp_max = (k_max << g_n) + 1;

	print("Testing n = ", n, ", p in [", zp_min, ", ", zp_max, "[.");

	print("For i = ", i_min, " to ", i_max);
	for (i = i_min, i_max - 1,

		\\ generate_primes
		for (id = 0, 2^log2_global_worksize - 1,
			j = (i << log2_global_worksize) + id;
			k = 15 * (j \ 8) + wheel[1 + j % 8];
			p = (k << g_n) + 1;
			if ((p % 3 == 0) || (p % 5 == 0), print("Error: wheel."); return;);
			if (prp(p),
		
				\\ init_factors
				a = 3;
				if (p % 3 != 2,
					a += 2;
					if (kro_vector[1 + 256 * ((5 - 3) / 2) + (p % 5)] >= 0,
						a += 2;
						if (kro_vector[1 + 256 * ((7 - 3) / 2) + (p % 7)] >= 0,
							a += 4;
							while (a < 256,
								if (kro_vector[1 + 256 * ((a - 3) / 2) + (p % a)] < 0, break;);
								a += 2;
							);
							if (a >= 256, a = 0;);
						);
					);
				);
					
				if (a != 0,
					\\ a^{(p - 1)/2} = -1 <=> a^{k*2^n} = -1. (a^k)^{2*j + 1} are the roots of b^{2^n} + 1 = 0 (mod p)
					c = Mod(a, p)^k;
					c0sq = c * c;
					
					\\ check_factors
					for (j = 1, 2^(n - 1),
						b = lift(c); if (b % 2 != 0, b = p - b;); \\ b is even
						if (b <= 2000000000,
							if (Mod(b, p)^(2^n) + 1 == 0,
								print(p, " | ", b, "^{2^", n, "}+1"); write("sieve.txt", p, " | ", b, "^{2^", n, "}+1");,
								if (isprime(p), print(p, " doesn't divide ", b, "^{2^", n, "}+1"); write("sieve.txt", p, " doesn't divide ", b, "^{2^", n, "}+1"););
							);
						);

						c *= c0sq;
					);
				);
			);
		);
	);
}

main()
