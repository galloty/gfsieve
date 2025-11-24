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
	kro_vector = vectorsmall(128 * 256);
	forstep (j = 3, 255, 2,
		for (i = 0, j - 1,
			kro_vector[1 + 256 * ((j - 3) / 2) + i] = kronecker(i, j);
		);
	);

	n = 23; g_n = n + 1;
	
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
	
	log2_global_worksize = 22; global_worksize = 2^log2_global_worksize;
	factors_loop = 2^10; N_2_factors_loop = 2^(n - 1) / factors_loop;

	p_min = 1; p_max = 2; unit = 1e15;

	f = unit * 8.0 / 15 / 2^(n + 1 + log2_global_worksize);
	i_min = floor(p_min * f); i_max = ceil(p_max * f);

	print("For i = ", i_min, " to ", i_max);
	for (i = i_min, i_max - 1,

		k_vector = List();

		\\ generate_primes
		for (id = 0, global_worksize - 1,
			j = (i << log2_global_worksize) + id;
			k = 15 * (j \ 8) + wheel[1 + j % 8];
			p = (k << g_n) + 1;
			if ((p % 3 == 0) || (p % 5 == 0), print("Error: wheel."); return;);
			if (prp(p), listput(~k_vector, k););
		);
		
		c_vector = vector(#k_vector);
		a2k_vector = vector(#k_vector);

		\\ init_factors
		for (id = 1, #k_vector,
			k = k_vector[id]; p = (k << g_n) + 1;

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
			
			c = Mod(a, p)^k;
			c_vector[id] = lift(c);
			a2k_vector[id] = lift(c^2);
		);

		\\ check_factors
		for (j = 0, N_2_factors_loop - 1,
			for (id = 1, #k_vector,

				if (c_vector[id] > 0,
					k = k_vector[id]; p = (k << g_n) + 1;
					c = Mod(c_vector[id], p);
					a2k = a2k_vector[id];
					
					for (l = 0, factors_loop - 1,
						b = lift(c); if (b % 2 != 0, b = p - b;);
						if (b <= 2000000000,
							if (Mod(b, p)^(2^n) + 1 == 0,
								print(p, " | ", b, "^{2^", n, "}+1"); write("sieve_gp.txt", p, " | ", b, "^{2^", n, "}+1");,
								if (isprime(p), print(p, " doesn't divide ", b, "^{2^", n, "}+1"); write("sieve_gp.txt", p, " doesn't divide ", b, "^{2^", n, "}+1"););
							);
						);

						c *= a2k;		\\ c = a^{(2*i + 1).k}
					);
					
					c_vector[id] = lift(c);
				);
			);
		);
	);
}

main()
