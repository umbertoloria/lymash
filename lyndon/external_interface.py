#!/bin/env python3
from sys import argv
from factorization import cfl, icfl, cfl_icfl, d_cfl, d_icfl, d_cfl_icfl

algs = {
	'cfl': cfl,
	'icfl': icfl,
	'cfl_icfl': cfl_icfl,
	'cfl_comb': d_cfl,
	'icfl_comb': d_icfl,
	'cfl_icfl_comb': d_cfl_icfl
}

alg = algs[argv[1]]

facts = alg(argv[2])
print(' '.join(facts))