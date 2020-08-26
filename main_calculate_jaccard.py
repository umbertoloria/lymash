from format import read_from_fasta
from jaccard import jaccard
from sliding_windows import ksliding


def explore_kmer_size():
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


def fixed_kmer_size():
	file1 = input("Tell me the first file> ")
	file2 = input("Tell me the second file> ")
	kmer_size = int(input("Tell me the kmer size (usually 21)> "))
	title1, read1 = read_from_fasta(file1)
	title2, read2 = read_from_fasta(file2)
	print("Comparing:", file1, "with", file2)
	print(file1, "'s title:", title1)
	print(file2, "'s title:", title2)
	print("kmer size:", kmer_size)
	print()

	p = jaccard(set(ksliding(read1, kmer_size)), set(ksliding(read2, kmer_size)))
	print("Similarity: %10.9f" % p)
