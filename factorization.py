from lyndon.factorization import cfl, icfl, cfl_icfl, d_cfl, d_icfl, d_cfl_icfl


def get_factors(alg, text):
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


def subdivide(text, cut_length):
	i = 0
	divisions = []
	while i + cut_length < len(text):
		division = text[i:i + cut_length]
		divisions.append(division)
		i += cut_length
	divisions.append(text[i:])
	return divisions


def get_fingers_after_split(alg, text, cut_length):
	fingers = []
	text_truncated = subdivide(text, cut_length)
	for text_part in text_truncated:
		factors = get_factors(alg, text_part)
		partial_fingers = [len(f) for f in factors]
		fingers.extend(partial_fingers)
	return tuple(fingers)
