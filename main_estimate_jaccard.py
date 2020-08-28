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
	for split in range(50, 301, 25):  # 50-300
		for window_size in range(3, 9):  # 3-8
			configs.append((split, window_size))

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
			"cfl_icfl": [],
			"cfl_comb": [],
			"icfl_comb": [],
			"cfl_icfl_comb": []
		}
	}

	for factorization_method in data["estimations"].keys():
		for split, window_size in configs:
			jaccard_estimated = estimation(factorization_method, read1, read2, split, window_size)
			diff = jaccard_calculated - jaccard_estimated

			if abs(diff) <= max_distance_to_show:
				coppie_famose[(split, window_size, factorization_method)] += 1
				data["estimations"][factorization_method].append({
					"split": split,
					"window size": window_size,
					"estimation": jaccard_estimated,
					"accuracy": diff
				})

	for output_function in output_functions:
		output_function(data)


def export_stdout(data):
	print("Deterministicly calculated (kmer size %3d): %10.9f" % (
		data["calculation"]["kmer size"], data["calculation"]["result"]))
	closest = 1
	for d in data["estimations"]["cfl"]:
		diff = d['accuracy']
		if abs(closest) > abs(diff):
			closest = abs(diff)
			print("( %3d - %2d ) %10.9f   (%5.4f) ***" % (d['split'], d['window size'], d['estimation'], d['accuracy']))
		else:
			print("( %3d - %2d ) %10.9f   (%5.4f)" % (d['split'], d['window size'], d['estimation'], d['accuracy']))


def export_csv(data):
	fixed = "___ split,___ window size,___ estimate,___ accuracy,"
	algs = data["estimations"].keys()
	with open(data["file1"] + '-' + data["file2"] + '.csv', mode='w') as file:
		file.write("kmer size,result,")
		for alg in algs:
			file.write(fixed.replace("___", alg))
		file.write("\n")

		i = {}
		for alg in algs:
			i[alg] = 1
		lens = {}
		for alg in algs:
			lens[alg] = len(data["estimations"][alg])

		file.write(str(data["calculation"]["kmer size"]) + "," + str(data["calculation"]["result"]) + ",")
		for alg in algs:
			if lens[alg] > 0:
				file.write(str(data["estimations"][alg][0]["split"]) + ",")
				file.write(str(data["estimations"][alg][0]["window size"]) + ",")
				file.write(str(data["estimations"][alg][0]["estimation"]) + ",")
				file.write(str(data["estimations"][alg][0]["accuracy"]) + ",")
			else:
				file.write(",,,,")
		file.write("\n")

		while True:
			stop = True
			for alg in algs:
				if i[alg] < lens[alg]:
					stop = False
					break
			if stop:
				break
			file.write(",,")
			for alg in algs:
				if i[alg] < lens[alg]:
					file.write(str(data["estimations"][alg][i[alg]]["split"]) + ",")
					file.write(str(data["estimations"][alg][i[alg]]["window size"]) + ",")
					file.write(str(data["estimations"][alg][i[alg]]["estimation"]) + ",")
					file.write(str(data["estimations"][alg][i[alg]]["accuracy"]) + ",")
				else:
					file.write(",,,,")
			file.write("\n")
			for alg in algs:
				i[alg] += 1
