import os
from factorization import get_factors, get_fingers_after_split, subdivide
from lyndon.utils import _complement
from sequences.Sequence import FastaSequence


def mash_interfacing_suite_main(action):
	if action == 5:
		preprocess_directory()
	elif action == 6:
		encode_fingerprints_in_ascii()


def preprocess_directory():
	dir = input("Tell me the directory with the file to preprocess> ")
	factorization = "cfl_icfl_comb"
	# input("Tell me the algorithm for factoring\n   (cfl,icfl,cfl_icfl,cfl_comb,icfl_comb,cfl_icfl_comb)> ")
	alf = " 0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZàèìòù%^-+:_çé!?£$§°*"

	if not os.path.isdir(dir + "_preprocessed"):
		os.mkdir(dir + "_preprocessed")
	for root, dirs, filepaths in os.walk(dir):
		for filepath in filepaths:
			if filepath.endswith(".fasta"):
				seq = FastaSequence(dir + "/" + filepath)
				comp_read = _complement(seq.get_data())
				print("read : {} reversed : {}".format(seq.get_data(), comp_read))

				fingerprint = [len(f) for f in get_factors(factorization, seq.get_data())]
				comp_fingerprint = [len(f) for f in get_factors(factorization, comp_read)]

				file = open(os.path.join(dir + "_preprocessed", "preprocessed_" + filepath), "w")
				file.write("@" + seq.get_title() + " preprocessed encoding lyndon fact lengths with ascii code\n")

				print(filepath)
				print(fingerprint)
				for finger in fingerprint:
					file.write(alf[finger])

				file.write("#")
				print(comp_fingerprint)
				for finger in comp_fingerprint:
					file.write(alf[finger])

				'''for fingers_row in subdivide(fingers, 70):
					for finger in fingers_row:
						num = ord('0') - 1 + finger
						if num > ord('Z'):
							print(" *** using not printable character (" + str(num) + "):", chr(num))
						file.write(chr(num))
					file.write("\n")'''
				file.close()


def encode_fingerprints_in_ascii():
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

	seq1 = FastaSequence(file1)
	seq2 = FastaSequence(file2)

	print("Translating and comparing:", file1, "with", file2, "using Mash")
	print(file1, "'s title:", seq1.get_title())
	print(file2, "'s title:", seq2.get_title())

	a, b = get_codified_factors(seq1.get_data(), seq2.get_data(), subdivision)

	af = open(file1 + ".translated", "w")
	af.write(seq1.get_title())
	af.write("\n")
	for part in subdivide(a, 70):
		af.write(part)
		af.write("\n")
	af.close()

	bf = open(file2 + ".translated", "w")
	bf.write(seq2.get_title())
	bf.write("\n")
	for part in subdivide(b, 70):
		bf.write(part)
		bf.write("\n")
	bf.close()
