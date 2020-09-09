from format import read_from_fasta
from jaccard import jaccard
from sliding_windows import ksliding
from datetime import datetime


def main_explore_kmer_size():
	file1 = input("Tell me the first file> ")
	file2 = input("Tell me the second file> ")
	title1, read1 = read_from_fasta(file1)
	title2, read2 = read_from_fasta(file2)
	print("Comparing:", file1, "with", file2)
	print(file1, "'s title:", title1)
	print(file2, "'s title:", title2)
	print()

	p = jaccard(set(ksliding(read1, 1)), set(ksliding(read2, 1)))
	print("%3d => %10.9f" % (1, p))

	prev_zero = 0
	tolerance = 10
	kmer_size = 2
	while prev_zero < tolerance:
		j = jaccard(set(ksliding(read1, kmer_size)), set(ksliding(read2, kmer_size)))
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


def calc_similarity_with_kmers_from_fasta_pair(file1, file2, kmer_size, *output_functions):
	title1, read1 = read_from_fasta(file1)
	title2, read2 = read_from_fasta(file2)
	before = datetime.now()
	similarity = jaccard(set(ksliding(read1, kmer_size)), set(ksliding(read2, kmer_size)))
	after = datetime.now()
	duration = after - before
	for output_function in output_functions:
		output_function(file1, file2, title1, title2, kmer_size, similarity, duration)


def main_fixed_kmer_size():
	def stdout_output_function(file1, file2, title1, title2, kmer_size, similarity, duration):
		# print("Comparing:", file1, "with", file2)
		# print(file1, "'s title:", title1)
		# print(file2, "'s title:", title2)
		# print("kmer size:", kmer_size, "\n")
		# print("Similarity: %10.9f" % similarity)
		n1 = file1.split("/")[-1].split(".")[0]
		n2 = file2.split("/")[-1].split(".")[0]
		print("%2s <-> %2s  (k=%d) -> %6.2f" % (n1, n2, kmer_size, similarity * 100) + "%" + (
				"    (%d microseconds)" % duration.microseconds))

	enable_pair_computation = input("Compute comparations from pairs (Y/n)> ")
	enable_pair_computation = True if enable_pair_computation == "" else (enable_pair_computation == "Y")

	if enable_pair_computation:

		import os
		print("Tell me files to compare ")
		files = []
		path = input("> ")
		while path != "":
			if not os.path.exists(path):
				print("File o directory non esistente")
			elif os.path.isdir(path):
				new_files = os.listdir(path)
				new_files.sort()
				for new_file in new_files:
					new_file = path + "/" + new_file
					if new_file not in files:
						files.append(new_file)
			else:
				if path not in files:
					files.append(path)
			path = input("> ")

		kmer_size = input("Tell me the kmer size (usually 21)> ")
		kmer_size = 21 if kmer_size == "" else int(kmer_size)

		from itertools import combinations
		for file1, file2 in combinations(files, 2):
			calc_similarity_with_kmers_from_fasta_pair(file1, file2, kmer_size, stdout_output_function)

	else:

		file1 = input("Tell me the first file> ")
		file2 = input("Tell me the second file> ")
		kmer_size = input("Tell me the kmer size (usually 21)> ")
		kmer_size = 21 if kmer_size == "" else int(kmer_size)
		calc_similarity_with_kmers_from_fasta_pair(file1, file2, kmer_size, stdout_output_function)
