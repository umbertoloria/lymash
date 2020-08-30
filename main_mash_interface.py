from factorization import get_fingers_after_split, subdivide
from format import read_from_fasta


def main():
	file1 = input("Tell me the first FASTA filepath to compare> ")
	file2 = input("Tell me the second FASTA filepath to compare> ")
	factorization = input(
		"Tell me the algorithm for factoring\n (cfl,icfl,cfl_icfl,cfl_comb,icfl_comb,cfl_icfl_comb)> ")
	subdivision = int(input("Tell me the subdivision lenght (max 200 for ASCII sake)> "))

	def get_codified_factors(a, b, split):
		a_factors_lengths = get_fingers_after_split(factorization, a, split)
		b_factors_lengths = get_fingers_after_split(factorization, b, split)
		alf = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
		lengths = list(set(a_factors_lengths).union(set(b_factors_lengths)))
		lengths.sort()
		# print("lunghezze veramente usate:", len(lengths))
		# print(a_factors_lengths)
		a_translated_factors = ""
		for length in a_factors_lengths:
			a_translated_factors += alf[lengths.index(length)]
		# print(b_factors_lengths)
		b_translated_factors = ""
		for length in b_factors_lengths:
			b_translated_factors += alf[lengths.index(length)]
		return a_translated_factors, b_translated_factors

	title1, read1 = read_from_fasta(file1)
	title2, read2 = read_from_fasta(file2)

	print("Translating and comparing:", file1, "with", file2, "using Mash")
	print(file1, "'s title:", title1)
	print(file2, "'s title:", title2)

	a, b = get_codified_factors(read1, read2, subdivision)

	af = open(file1 + ".translated", "w")
	af.write(title1)
	af.write("\n")
	for part in subdivide(a, 70):
		af.write(part)
		af.write("\n")
	af.close()

	bf = open(file2 + ".translated", "w")
	bf.write(title1)
	bf.write("\n")
	for part in subdivide(b, 70):
		bf.write(part)
		bf.write("\n")
	bf.close()
