from collections import defaultdict

from lyndon.utils import _complement
from sequences.Sequence import Sequence
from traditional_technique import jaccard, ksliding, jaccard_on_kmers

FACTORIZATIONS = ['cfl', 'icfl', 'cfl_icfl', 'cfl_comb', 'icfl_comb', 'cfl_icfl_comb']


def factorize(factorization, text, split=0):
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
	factors = []
	if split == 0:
		factors = func(text)
	else:
		for text_part in subdivide(text, split):
			factors += func(text_part)
	return tuple(factors)


def jaccard_on_kfingers(str1: str, str2: str, factorization: str, k: int, use_super_fp: bool, split: int = 0):
	afingerprint = tuple(create_fingerprint(str1, factorization, use_super_fp, split))
	bfingerprint = tuple(create_fingerprint(str2, factorization, use_super_fp, split))
	akfingers = set(ksliding(afingerprint, k, False))
	bkfingers = set(ksliding(bfingerprint, k, False))
	return jaccard(akfingers, bkfingers)


# TODO: questa funzione non usa la super-fingerprint!
def new_combined_technique_analyzer(seq1: Sequence, seq2: Sequence, kmer_size: int, tolerance: float, *out_funcs):
	data = JaccardFactResult()
	data.set_sequence1(seq1)
	data.set_sequence2(seq2)

	jaccard_calculated = jaccard_on_kmers(seq1.get_data(), seq2.get_data(), kmer_size)
	data.add_calculated_jaccard(kmer_size, jaccard_calculated)

	result = {'jaccard': jaccard_calculated, 'factorizations': {}}

	for factorization in FACTORIZATIONS:
		sopravvissuti = []

		for split in range(50, 301, 25):

			a_normal_fingerprint = tuple([len(f) for f in factorize(factorization, seq1.get_data(), split)])
			b_normal_fingerprint = tuple([len(f) for f in factorize(factorization, seq2.get_data(), split)])

			for window_size in range(3, 9):  # 3-8
				if window_size > len(a_normal_fingerprint) or window_size > len(b_normal_fingerprint):
					continue
				a_kfingers = set(ksliding(a_normal_fingerprint, window_size))
				b_kfingers = set(ksliding(b_normal_fingerprint, window_size))
				jaccard_estimated = jaccard(a_kfingers, b_kfingers)

				# Instead of using the code below, we separate the fingerprint calculation from
				# the kfingers extrapolation. This way we don't keep re-calculating it for every window_size value.
				# jaccard_on_kfingers(seq1.get_data(), seq2.get_data(), factorization, window_size, False, split)

				diff = jaccard_calculated - jaccard_estimated
				if abs(diff) <= tolerance:
					sopravvissuti.append((split, window_size))
					data.add_estimated_jaccard(factorization, split, window_size, jaccard_estimated, diff)

		result['factorizations'][factorization] = sopravvissuti

	for out_func in out_funcs:
		out_func(data)

	return result


def create_fingerprint(text, factorization: str, use_super_fp: bool, split: int = 0):
	fingerprint = [len(f) for f in factorize(factorization, text, split)]
	if use_super_fp:
		fingerprint += [0]
		fingerprint += [len(f) for f in factorize(factorization, _complement(text), split)]
	return fingerprint


def subdivide(text, cut_length):
	i = 0
	divisions = []
	while i + cut_length < len(text):
		division = text[i:i + cut_length]
		divisions.append(division)
		i += cut_length
	divisions.append(text[i:])
	return divisions


class JaccardFactResult:

	def __init__(self):
		self.__calculations = {}
		self.__estimations = defaultdict(dict)

	def get_sequence1(self) -> Sequence:
		return self.__seq1

	def set_sequence1(self, seq: Sequence):
		self.__seq1: Sequence = seq

	def get_sequence2(self) -> Sequence:
		return self.__seq2

	def set_sequence2(self, seq: Sequence):
		self.__seq2: Sequence = seq

	def add_calculated_jaccard(self, kmer_size: int, jacc: float):
		self.__calculations[kmer_size] = jacc

	def get_calculated_jaccard(self, kmer_size: int) -> float:
		return self.__calculations[kmer_size]

	def add_estimated_jaccard(self, factorization: str, split: int, window_size: int, jacc: float, accuracy: float):
		self.__estimations[factorization][(split, window_size)] = {'estimation': jacc, 'accuracy': accuracy}

	def get_factorizations_used(self):
		return self.__estimations.keys()

	def get_estimated_jaccard_from_factorization(self, factorization: str):
		return self.__estimations[factorization]


def get_csv_exporter(dirname):
	def export_csv(data: JaccardFactResult):
		fixed = '___ split,___ window size,___ estimate,___ accuracy,'

		csvname = data.get_sequence1().get_name() + '-' + data.get_sequence2().get_name() + '.csv'

		kmer_size = 21

		estimations = {}
		for factorization in FACTORIZATIONS:
			estimations[factorization] = []
			ests = data.get_estimated_jaccard_from_factorization(factorization)
			for (split, window_size), value in ests.items():
				estimations[factorization].append(
					[str(split), str(window_size), str(value['estimation']), str(value['accuracy'])])

		with open(dirname + '/' + csvname, 'w') as csv:
			csv.write('kmer size,result,')
			for factorization in FACTORIZATIONS:
				csv.write(fixed.replace('___', factorization))
			csv.write('\n')

			i = {}
			for factorization in FACTORIZATIONS:
				i[factorization] = 1
			lens = {}
			for factorization in FACTORIZATIONS:
				lens[factorization] = len(estimations[factorization])

			csv.write(str(kmer_size) + ',' + str(data.get_calculated_jaccard(kmer_size)) + ',')
			for factorization in FACTORIZATIONS:
				if lens[factorization] > 0:
					csv.write(','.join(estimations[factorization][0]) + ',')
				else:
					csv.write(',,,,')
			csv.write('\n')

			while True:
				stop = True
				for factorization in FACTORIZATIONS:
					if i[factorization] < lens[factorization]:
						stop = False
						break
				if stop:
					break
				csv.write(',,')
				for factorization in FACTORIZATIONS:
					if i[factorization] < lens[factorization]:
						csv.write(",".join(estimations[factorization][i[factorization]]) + ",")
					else:
						csv.write(",,,,")
				csv.write("\n")
				for factorization in FACTORIZATIONS:
					i[factorization] += 1

	return export_csv


def get_grafico_exporter(data: JaccardFactResult):
	import matplotlib.pyplot as plt
	plt.figure(num=None, figsize=(15, 10), dpi=80, facecolor="w", edgecolor="k")
	values = []
	for factorization in FACTORIZATIONS:
		val = 0
		ests = data.get_estimated_jaccard_from_factorization(factorization)
		for (split, window_size), value in ests.items():
			val += 1
		values.append(val)
	plt.bar(FACTORIZATIONS, values)
	# plt.axis([0, 1, 0, 1])

	kmer_size = 21
	graphname = data.get_sequence1().get_name() + "-" + data.get_sequence2().get_name() + ".csv"

	plt.xlabel("Fattorizzazioni")
	plt.ylabel("Precisione (oracolo: " + str(data.get_calculated_jaccard(kmer_size)) + ")")
	plt.title(graphname)
	plt.grid()
	plt.show()
