def get_factors(factorization, text):
	from lyndon.factorization import cfl, icfl, cfl_icfl, d_cfl, d_icfl, d_cfl_icfl
	factorizations = {
		'cfl': cfl,
		'icfl': icfl,
		'cfl_icfl': cfl_icfl,
		'cfl_comb': d_cfl,
		'icfl_comb': d_icfl,
		'cfl_icfl_comb': d_cfl_icfl
	}

	func = factorizations[factorization]
	return tuple(func(text))

# TODO: use cfl_max...
# TODO: MERGE WITH NEW_COMBINED_TECHNIQUE.py

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
