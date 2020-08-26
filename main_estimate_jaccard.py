from factorization import get_fingers_after_split
from format import read_from_fasta
from jaccard import jaccard
from sliding_windows import ksliding


def estimate_jaccard_from_filepaths(file1, file2, max_distance_to_show, coppie_famose, *output_functions):
	def estimation(factorization_method, a, b, split, window_size):
		a_factors_lengths = get_fingers_after_split(factorization_method, a, split)
		a_fingerprint = set(ksliding(a_factors_lengths, window_size))
		b_factors_lengths = get_fingers_after_split(factorization_method, b, split)
		b_fingerprint = set(ksliding(b_factors_lengths, window_size))
		return jaccard(a_fingerprint, b_fingerprint)

	title1, read1 = read_from_fasta(file1)
	title2, read2 = read_from_fasta(file2)

	# print("Comparing:", file1, "with", file2)
	# print(file1, "'s title:", title1)
	# print(file2, "'s title:", title2)
	# print()

	jaccard_calculated = jaccard(set(ksliding(read1, 21)), set(ksliding(read2, 21)))

	configs = []
	for split in range(10, 301, 10):
		for window_size in range(2, 9):
			configs.append((split, window_size))

	print()
	data = {
		"file1": file1.split("/")[-1],
		"file2": file2.split("/")[-1],
		"calculation": {
			"kmer size": 21,
			"result": jaccard_calculated
		},
		"estimations": {
			"cfl": [],
			"icfl": [],
		}
	}

	for factorization_method in ["cfl", "icfl"]:
		for split, window_size in configs:
			jaccard_estimated = estimation(factorization_method, read1, read2, split, window_size)
			diff = jaccard_calculated - jaccard_estimated

			if abs(diff) <= max_distance_to_show:
				coppie_famose[(split, window_size)] += 1
				data["estimations"][factorization_method].append({
					"split": split,
					"window size": window_size,
					"result": jaccard_estimated,
					"accuracy": diff
				})

	for output_function in output_functions:
		output_function(data)


'''
def main():
	file1 = input("Tell me the first FASTA filepath to compare> ")
	file2 = input("Tell me the second FASTA filepath to compare> ")
	# factorization_method = input(
	#		"Tell me the algorithm for factoring\n   (cfl,icfl,cfl_icfl,cfl_comb,icfl_comb,cfl_icfl_comb)> ")
	if file1 == "":
		file1 = "test/1.fasta"
	if file2 == "":
		file2 = "test/2.fasta"
	# if factorization_method == "":
	#		factorization_method = "cfl"
	estimate_jaccard_from_filepaths(file1, file2, 0.05, ..., export_csv)
'''


def export_stdout(data):
	print("Deterministicly calculated (kmer size %3d): %10.9f" % (
		data["calculation"]["kmer size"], data["calculation"]["result"]))
	closest = 1
	for d in data["estimations"]["cfl"]:
		diff = d['accuracy']
		if abs(closest) > abs(diff):
			closest = abs(diff)
			print("( %3d - %2d ) %10.9f   (%5.4f) ***" % (d['split'], d['window size'], d['result'], d['accuracy']))
		else:
			print("( %3d - %2d ) %10.9f   (%5.4f)" % (d['split'], d['window size'], d['result'], d['accuracy']))


def export_csv(data):
	import csv

	fieldnames = ['kmer size', 'result',
	              'cfl split', 'cfl window size', 'cfl result', 'cfl accuracy',
	              'icfl split', 'icfl window size', 'icfl result', 'icfl accuracy', ]
	with open(data["file1"] + '-' + data["file2"] + '.csv', mode='w') as file:
		writer = csv.DictWriter(file, fieldnames=fieldnames)
		writer.writeheader()
		writer.writerow({
			'kmer size': data["calculation"]["kmer size"],
			'result': data["calculation"]["result"],
			'cfl split': data["estimations"]["cfl"][0]["split"]
			if len(data["estimations"]["cfl"]) > 0 else '',
			'cfl window size': data["estimations"]["cfl"][0]["window size"]
			if len(data["estimations"]["cfl"]) > 0 else '',
			'cfl result': data["estimations"]["cfl"][0]["result"]
			if len(data["estimations"]["cfl"]) > 0 else '',
			'cfl accuracy': data["estimations"]["cfl"][0]["accuracy"]
			if len(data["estimations"]["cfl"]) > 0 else '',
			'icfl split': data["estimations"]["icfl"][0]["split"]
			if len(data["estimations"]["icfl"]) > 0 else '',
			'icfl window size': data["estimations"]["icfl"][0]["window size"]
			if len(data["estimations"]["icfl"]) > 0 else '',
			'icfl result': data["estimations"]["icfl"][0]["result"]
			if len(data["estimations"]["icfl"]) > 0 else '',
			'icfl accuracy': data["estimations"]["icfl"][0]["accuracy"]
			if len(data["estimations"]["icfl"]) > 0 else '',
		})

		cfl_i = 1
		icfl_i = 1
		cfl_len = len(data["estimations"]["cfl"])
		icfl_len = len(data["estimations"]["icfl"])

		while cfl_i < cfl_len and icfl_i < icfl_len:
			writer.writerow({
				'kmer size': '',
				'result': '',
				'cfl split': data["estimations"]["cfl"][cfl_i]["split"],
				'cfl window size': data["estimations"]["cfl"][cfl_i]["window size"],
				'cfl result': data["estimations"]["cfl"][cfl_i]["result"],
				'cfl accuracy': data["estimations"]["cfl"][cfl_i]["accuracy"],
				'icfl split': data["estimations"]["icfl"][icfl_i]["split"],
				'icfl window size': data["estimations"]["icfl"][icfl_i]["window size"],
				'icfl result': data["estimations"]["icfl"][icfl_i]["result"],
				'icfl accuracy': data["estimations"]["icfl"][icfl_i]["accuracy"],
			})
			cfl_i += 1
			icfl_i += 1
		if cfl_i < cfl_len:
			writer.writerow({
				'kmer size': '',
				'result': '',
				'cfl split': data["estimations"]["cfl"][cfl_i]["split"],
				'cfl window size': data["estimations"]["cfl"][cfl_i]["window size"],
				'cfl result': data["estimations"]["cfl"][cfl_i]["result"],
				'cfl accuracy': data["estimations"]["cfl"][cfl_i]["accuracy"],
				'icfl split': '',
				'icfl window size': '',
				'icfl result': '',
				'icfl accuracy': '',
			})
			cfl_i += 1
		elif icfl_i < icfl_len:
			writer.writerow({
				'kmer size': '',
				'result': '',
				'cfl split': '',
				'cfl window size': '',
				'cfl result': '',
				'cfl accuracy': '',
				'icfl split': data["estimations"]["icfl"][icfl_i]["split"],
				'icfl window size': data["estimations"]["icfl"][icfl_i]["window size"],
				'icfl result': data["estimations"]["icfl"][icfl_i]["result"],
				'icfl accuracy': data["estimations"]["icfl"][icfl_i]["accuracy"],
			})
			icfl_i += 1
