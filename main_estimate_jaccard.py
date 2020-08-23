from factorization import get_fingers_after_text_subdividing
from format import read_from_fasta
from jaccard import jaccard
from sliding_windows import ksliding


def main():
	file1 = input("Tell me the first FASTA filepath to compare> ")
	file2 = input("Tell me the second FASTA filepath to compare> ")
	factorization_method = input(
		"Tell me the algorithm for factoring\n   (cfl,icfl,cfl_icfl,cfl_comb,icfl_comb,cfl_icfl_comb)> ")

	def estimation(a, b, subdivision, slide_length):
		a_factors_lengths = get_fingers_after_text_subdividing(factorization_method, a, subdivision)
		a_fingerprint = set(ksliding(a_factors_lengths, slide_length))
		b_factors_lengths = get_fingers_after_text_subdividing(factorization_method, b, subdivision)
		b_fingerprint = set(ksliding(b_factors_lengths, slide_length))
		return jaccard(a_fingerprint, b_fingerprint)

	title1, read1 = read_from_fasta(file1)
	title2, read2 = read_from_fasta(file2)

	print("Comparing:", file1, "with", file2)
	print(file1, "'s title:", title1)
	print(file2, "'s title:", title2)
	print()

	j = jaccard(set(ksliding(read1, 21)), set(ksliding(read2, 21)))

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
