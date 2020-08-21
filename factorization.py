import subprocess
from lyndon.factorization import cfl, icfl, cfl_icfl, d_cfl, d_icfl, d_cfl_icfl

def get_factors(alg, text):
	#cmd = ["../lyndon-master/python/external_interface.py", alg, text]
	#result = subprocess.run(cmd, stdout=subprocess.PIPE)
	algs = {
		'cfl': cfl,
		'icfl': icfl,
		'cfl_icfl': cfl_icfl,
		'cfl_comb': d_cfl,
		'icfl_comb': d_icfl,
		'cfl_icfl_comb': d_cfl_icfl
	}

	fun = algs[alg]
	return tuple(fun(text))
	#out = result.stdout.decode("utf-8")
	#if out.endswith('\n'):
	#	out = out[:-1]
	#return tuple(out.split(" "))


def subdivide(text, max_length):
	i = 0
	divisions = []
	while i + max_length < len(text):
		division = text[i:i + max_length]
		divisions.append(division)
		i += max_length
	divisions.append(text[i:])
	return divisions


def get_fingers_after_text_subdividing(alg, text, max_length):
	fingers = []
	text_truncated = subdivide(text, max_length)
	for text_part in text_truncated:
		factors = get_factors(alg, text_part)
		partial_fingers = [len(f) for f in factors]
		fingers.extend(partial_fingers)
	return tuple(fingers)
