#!/bin/env python3

from factorization import get_fingers_after_text_subdividing, subdivide
from format import get_reads_from_fq_file, read_from_fasta
from jaccard import jaccard
from sliding_windows import ksliding


def not_estimation(a, b, kmer_size):
	a_kmers = ksliding(a, kmer_size)
	b_kmers = ksliding(b, kmer_size)

	j = jaccard(set(a_kmers), set(b_kmers))

	# print("KMER SIZE:", kmer_size)
	# print("JACCARD SIMILARITY: %10.9f" % j)
	# print()

	return j


def estimation(a, b, truncate, slide_length):
	a_factors_lengths = tuple(get_fingers_after_text_subdividing("cfl", a, truncate))
	b_factors_lengths = tuple(get_fingers_after_text_subdividing("cfl", b, truncate))

	a_fingerprint = ksliding(a_factors_lengths, slide_length)
	b_fingerprint = ksliding(b_factors_lengths, slide_length)

	j = jaccard(set(a_fingerprint), set(b_fingerprint))

	# print("TRUNCATE:", truncate)
	# print("SLIDE LENGTH:", slide_length)
	# print("JACCARD SIMILARITY OF FINGERPRINTS: %10.9f", j)
	# print()

	return j


def from_fq_to_multiple_fastas(fq_filepath):
	reads = get_reads_from_fq_file(fq_filepath)

	i = 1
	for title, read in reads:
		f = open("in/" + str(i) + ".fasta", "w")
		i += 1
		parts = subdivide(read, 70)
		f.write(">")
		f.write(title)
		for part in parts:
			f.write("\n")
			f.write(part)
		f.close()


def main():
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

	'''
	j = not_estimation(s, t, 21)
	
	# truncate = int(argv[1])
	# slide_length = int(argv[2])

	configs = []
	for truncate in range(10, 301, 10):
		for slide_length in range(2, 15):
			configs.append(
				(truncate, slide_length)
			)

	closest = 1

	values = {}

	for truncate, slide_length in configs:
		est = estimation(s, t, truncate, slide_length)
		values[(truncate, slide_length)] = est
		diff = abs(est - j)
		if closest > diff:
			closest = diff
			print("( %3d - %2d ) %10.9f   (%5.4f)" % (truncate, slide_length, est, diff))
		else:
			print("( %3d - %2d ) %10.9f" % (truncate, slide_length, est))
	'''


if __name__ == "__main__":
	main()
