from collections import defaultdict

from Sequence import Sequence
from factorization import get_fingers_after_split
from jaccard import jaccard
from sliding_windows import ksliding

FACTORIZATIONS = ["cfl", "icfl", "cfl_icfl", "cfl_comb", "icfl_comb", "cfl_icfl_comb"]


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
		self.__estimations[factorization][(split, window_size)] = {
			"estimation": jacc,
			"accuracy": accuracy
		}

	def get_factorizations_used(self):
		return self.__estimations.keys()

	def get_estimated_jaccard_from_factorization(self, factorization: str):
		return self.__estimations[factorization]


def jaccard_thanks_factorizations(seq1: Sequence, seq2: Sequence, kmer_size: int, tolerance: float, *out_funcs):
	read1 = seq1.get_data()
	read2 = seq2.get_data()

	data = JaccardFactResult()
	data.set_sequence1(seq1)
	data.set_sequence2(seq2)

	jaccard_calculated = jaccard(set(ksliding(read1, kmer_size)), set(ksliding(read2, kmer_size)))
	data.add_calculated_jaccard(kmer_size, jaccard_calculated)

	result = {
		"jaccard": jaccard_calculated,
		"factorizations": {}
	}

	for factorization in FACTORIZATIONS:
		sopravvissuti = []

		for split in range(50, 301, 25):  # 50-300

			a_factors_lengths = get_fingers_after_split(factorization, read1, split)
			b_factors_lengths = get_fingers_after_split(factorization, read2, split)

			for window_size in range(3, 9):  # 3-8

				a_fingerprint = set(ksliding(a_factors_lengths, window_size))
				b_fingerprint = set(ksliding(b_factors_lengths, window_size))
				jaccard_estimated = jaccard(a_fingerprint, b_fingerprint)

				diff = jaccard_calculated - jaccard_estimated
				if abs(diff) <= tolerance:
					sopravvissuti.append((split, window_size))
					data.add_estimated_jaccard(factorization, split, window_size, jaccard_estimated, diff)

		result["factorizations"][factorization] = sopravvissuti

	for out_func in out_funcs:
		out_func(data)

	return result


'''
FIXME
def export_stdout(data: JaccardFactResult):
	kmer_size = 21

	print("Deterministicly calculated (kmer size %3d): %10.9f" % (kmer_size, data.get_calculated_jaccard(kmer_size)))
	closest = 1
	for d in data.get_estimated_jaccard_from_factorization("cfl"):
		diff = d['accuracy']
		if abs(closest) > abs(diff):
			closest = abs(diff)
			print("( %3d - %2d ) %10.9f   (%5.4f) ***" % (d['split'], d['window size'], d['estimation'], d['accuracy']))
		else:
			print("( %3d - %2d ) %10.9f   (%5.4f)" % (d['split'], d['window size'], d['estimation'], d['accuracy']))
'''


def get_csv_exporter(dirname):
	def export_csv(data: JaccardFactResult):
		fixed = "___ split,___ window size,___ estimate,___ accuracy,"

		csvname = data.get_sequence1().get_name() + "-" + data.get_sequence2().get_name() + ".csv"

		kmer_size = 21

		estimations = {}
		for factorization in FACTORIZATIONS:
			estimations[factorization] = []
			ests = data.get_estimated_jaccard_from_factorization(factorization)
			for (split, window_size), value in ests.items():
				estimations[factorization].append(
					[str(split), str(window_size), str(value["estimation"]), str(value["accuracy"])])

		with open(dirname + "/" + csvname, 'w') as csv:
			csv.write("kmer size,result,")
			for factorization in FACTORIZATIONS:
				csv.write(fixed.replace("___", factorization))
			csv.write("\n")

			i = {}
			for factorization in FACTORIZATIONS:
				i[factorization] = 1
			lens = {}
			for factorization in FACTORIZATIONS:
				lens[factorization] = len(estimations[factorization])

			csv.write(str(kmer_size) + "," + str(data.get_calculated_jaccard(kmer_size)) + ",")
			for factorization in FACTORIZATIONS:
				if lens[factorization] > 0:
					csv.write(",".join(estimations[factorization][0]) + ",")
				else:
					csv.write(",,,,")
			csv.write("\n")

			while True:
				stop = True
				for factorization in FACTORIZATIONS:
					if i[factorization] < lens[factorization]:
						stop = False
						break
				if stop:
					break
				csv.write(",,")
				for factorization in FACTORIZATIONS:
					if i[factorization] < lens[factorization]:
						csv.write(",".join(estimations[factorization][i[factorization]]) + ",")
					else:
						csv.write(",,,,")
				csv.write("\n")
				for factorization in FACTORIZATIONS:
					i[factorization] += 1

	return export_csv
