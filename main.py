#!/bin/env python3

from factorization import get_fingers_after_text_subdividing, subdivide
from format import get_reads_from_fq_file, read_from_fasta
from jaccard import jaccard
from sliding_windows import ksliding
from preprocess import preprocess


def from_fq_to_multiple_fastas(fq_filepath):
	reads = get_reads_from_fq_file(fq_filepath)
	i = 1
	for title, read in reads:
		print("creating file " + str(i))
		f = open("in/" + str(i) + ".fasta", "w")
		i += 1
		parts = subdivide(read, 70)
		f.write(">")
		f.write(title)
		for part in parts:
			f.write("\n")
			f.write(part)
		f.close()


def not_estimation(a, b, kmer_size):
	a_kmers = ksliding(a, kmer_size)
	b_kmers = ksliding(b, kmer_size)
	j = jaccard(set(a_kmers), set(b_kmers))
	# print("KMER SIZE:", kmer_size)
	# print("JACCARD SIMILARITY: %10.9f" % j)
	# print()
	return j


def main_calculate_jaccard():
	computation = {
		"pairs": (
			("in/1.fasta", "in/2.fasta"),
			("in/1.fasta", "in/3.fasta"),
			("in/1.fasta", "in/4.fasta"),
			("in/1.fasta", "in/5.fasta"),
		)
	}

	for file1, file2 in computation["pairs"]:

		title1, read1 = read_from_fasta(file1)
		title2, read2 = read_from_fasta(file2)

		print("Comparing:", file1, "with", file2)
		print(file1, "'s title:", title1)
		print(file2, "'s title:", title2)
		print()

		j = not_estimation(read1, read2, 1)
		p = j
		print("%3d => %10.9f" % (1, j))

		prev_zero = 0
		tolerance = 10
		kmer_size = 2
		while prev_zero < tolerance:
			j = not_estimation(read1, read2, kmer_size)
			if p > 0:
				diff = 100 - j / p * 100
				prev_zero = 0
			else:
				diff = 0
				prev_zero += 1
			p = j
			print("%3d => %10.9f   (diminuisce del %5.3f" % (kmer_size, j, diff) + '%)')
			kmer_size += 1

		print()
		print()


def main_estimate_jaccard(file1, file2, factorization_method):
	def estimation(a, b, subdivision, slide_length):
		print("1 read: " + factorization_method + " + fingerprint")
		a_factors_lengths = get_fingers_after_text_subdividing(factorization_method, a, subdivision)
		a_fingerprint = set(ksliding(a_factors_lengths, slide_length))
		print("2 read: " + factorization_method + " + fingerprint")
		b_factors_lengths = get_fingers_after_text_subdividing(factorization_method, b, subdivision)
		b_fingerprint = set(ksliding(b_factors_lengths, slide_length))
		print("Calculating Jaccard similarity from fingerprints")
		j = jaccard(a_fingerprint, b_fingerprint)
		print("SUBDIVISION:", subdivision)
		print("SLIDE LENGTH:", slide_length)
		print("JACCARD SIMILARITY OF FINGERPRINTS: %10.9f", j)
		print()
		return j

	title1, read1 = read_from_fasta(file1)
	title2, read2 = read_from_fasta(file2)

	print("Comparing:", file1, "with", file2)
	print(file1, "'s title:", title1)
	print(file2, "'s title:", title2)
	print()

	j = not_estimation(read1, read2, 21)

	configs = []
	for truncate in range(10, 301, 10):
		for slide_length in range(2, 15):
			configs.append((truncate, slide_length))

	# Starting from subdivision=10, it is going to make way too much external terminal calls.
	# TODO: Internalize the Lyndon factorizer.
	closest = 1
	values = {}
	for truncate, slide_length in configs:
		est = estimation(read1, read2, truncate, slide_length)
		values[(truncate, slide_length)] = est
		diff = abs(est - j)
		if closest > diff:
			closest = diff
			print("( %3d - %2d ) %10.9f   (%5.4f)" % (truncate, slide_length, est, diff))
		else:
			print("( %3d - %2d ) %10.9f" % (truncate, slide_length, est))


if __name__ == "__main__":
	print("1: Export FASTA from FQ")
	print("2: Calculate Jaccard similarity")
	print("3: Estimate Jaccard similarity")
	print("4: Preprocess all the files in a directory")
	print("5: Quit")
	print()
	action = int(input("Tell me what to do > "))

	if action == 1:
		filepath = input("Tell me the filepath of the FQ file to \"explode\" > ")
		from_fq_to_multiple_fastas(filepath)
	elif action == 2:
		main_calculate_jaccard()
	elif action == 3:
		file1 = input("Tell me the first FASTA filepath to compare> ")
		file2 = input("Tell me the second FASTA filepath to compare> ")
		algo = input("Tell me the algorithm for factoring\n (cfl,icfl,cfl_icfl,cfl_comb,icfl_comb,cfl_icfl_comb)> ")
		main_estimate_jaccard(file1, file2, algo)
	elif action == 4:
		dir = input("Tell me the directory with the file to preprocess>")
		algo = input("Tell me the algorithm for factoring\n (cfl,icfl,cfl_icfl,cfl_comb,icfl_comb,cfl_icfl_comb)> ")
		subdivision = input("Tell me the subdivision lenght (max 200 for ASCII sake)")
		preprocess(dir, algo, int(subdivision))
	else:
		print("Goodbye...️")
