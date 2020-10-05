#!/bin/env python3

import main_export_fq_in_fasta_files
from files_manager import input_files
from traditional_technique import jaccard_on_kmers
from suites import jaccard_kmers_suite, mash_interfacing_suite, jaccard_kfingers_suite

if __name__ == "__main__":
	print("FORMATS")
	print("1: Export FASTA from FQ")
	print("")
	print("CALCULATE SIMILARITY: JACCARD ON K-MERS")
	print("2: Calculate Jaccard similarity (exploring various kmer size)")
	print("3: Calculate Jaccard similarity (fixed kmer size)")
	print("")
	print("ESTIMATE SIMILARITY: JACCARD ON K-FINGERS")
	print("4: Estimate Jaccard similarity (estimate on pairs of some files)")
	print("")
	print("MASH INTERFACING")
	print("5: Preprocess all the files in a directory")
	print("6: Encode in ASCII the fingers of two FASTA files")
	print("7: WE (pure a me plz)")
	print("8: Antonio dammi un nome...")
	print("9: Crea dataset di long long reads")
	print("10: Crea dataset generico")
	print("11: Dataset short reads scrivi, kmer, kfinger")
	print("")
	print("BOH...")
	print("_: Quit")
	print("")
	action = int(input("\nTell me what to do > "))

	# TODO: remove numbers for suites
	if 2 <= action <= 3:
		jaccard_kmers_suite.jaccard_kmers_suite_main(action)
	elif action == 4:
		jaccard_kfingers_suite.jaccard_kfingers_suite_main(action)
	elif 5 <= action <= 6:
		mash_interfacing_suite.mash_interfacing_suite_main(action)

	if action == 1:
		main_export_fq_in_fasta_files.main()

	elif action == 7:
		import os
		from sequences.Sequence import FastaSequence
		import matplotlib.pyplot as plt

		ascisse = []
		ascisse2 = []
		ascisse3 = []
		ordinate = []
		ordinate2 = []
		ordinate3 = []
		from main_estimate_jaccard import estimate_jaccard

		# CLI inputs
		files = input_files()
		kmer_size = input("Tell me the kmer size (21) > ")
		kmer_size = 21 if kmer_size == "" else int(kmer_size)
		print()
		factorization = input("Tell me the factorization (cfl_icfl_comb)")
		if factorization == "":
			factorization = "cfl_icfl_comb"
		kfinger_size = input("tell me the kfinger sizes(3)")
		kfinger_size = 3 if kfinger_size == "" else int(kfinger_size)

		use_super_fp = input("Tell me if you want the fingerprint to have [...-0-rev_compls] [Y/n]> ")
		use_super_fp = (use_super_fp == "Y" or use_super_fp == "y" or use_super_fp == "")

		max = 0
		import copy

		for i in range(len(files)):
			file = files[i]
			seq1 = FastaSequence(file)
			read1 = seq1.get_data()
			from factorization import get_factors
			from lyndon.utils import _complement

			factors = [len(f) for f in get_factors(factorization, read1)]
			factors2 = [len(f) for f in get_factors(factorization, _complement(read1))]
			factors.sort()
			factors2.sort()
			print(factors[-1])
			print("max = {}".format(max))
			if factors[-1] > max:
				max = copy.deepcopy(factors[-1])
			if factors[-1] > max:
				max = copy.deepcopy(factors2[-1])

		print("il max è {}".format(max))
		y = 0

		files = sorted(files, key=lambda x: int(x.split(".")[0].split("/")[1]))

		for i in range(len(files) - 1):
			file1 = files[i]
			file2 = files[i + 1]

			print("Comparing", file1, "-", file2)
			seq1 = FastaSequence(file1)
			seq2 = FastaSequence(file2)

			calc = jaccard_on_kmers(seq1.get_data(), seq2.get_data(), kmer_size)
			estim = estimate_jaccard(seq1, seq2, factorization, kfinger_size, use_super_fp)
			print(calc, estim)

			import subprocess

			out = subprocess.Popen(["mash", "dist", file1, file2, "-s", "1000"], stdout=subprocess.PIPE,
			                       stderr=subprocess.STDOUT)
			stdout, stderr = out.communicate()
			output = str(stdout).split('t')[-1].split('\\')[0]
			output = output.split("/")
			mashresult1 = int(output[0]) / int(output[1])

			file1pre = file1.replace("/", "_preprocessed/preprocessed_")
			file2pre = file2.replace("/", "_preprocessed/preprocessed_")

			out = subprocess.Popen(["mash", "dist", file1pre, file2pre, "-s", "1000", "-k", "3", "-z",
			                        "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZàèìòù%^-+:_çé!?£$§°*"],
			                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			stdout, stderr = out.communicate()
			output = str(stdout).split('t')[-1].split('\\')[0]
			output = output.split("/")
			mashresult2 = int(output[0]) / int(output[1])

			# ordinate2.append(abs(mashresult1-mashresult2))
			# ascisse2.append(mashresult1)
			ordinate2.append(mashresult1)
			ascisse2.append(calc)

			ordinate3.append(mashresult2)
			ascisse3.append(calc)

			ordinate.append(estim)
			ascisse.append(calc)

			print("---------------------------------")
			print()
			if abs(calc - estim) > 0.25:
				y += 1

		# plt.plot(ascisse,	ordinate)
		plt.plot(ascisse2, ordinate2, "b.")
		plt.plot(ascisse, ordinate, "r+")

		plt.plot(ascisse3, ordinate3, "y*")
		plt.axis([0, 1, 0, 1])
		plt.xlabel("Similarità stringhe")
		plt.ylabel("Accuracy^-1")
		plt.legend(["Mash k-mer", "Jaccard k-finger", "Mash k-finger"])
		plt.title("CIAO")
		plt.grid()
		plt.show()

		print(y)

	elif action == 8:
		import os
		from collections import defaultdict
		from itertools import combinations
		from sequences.Sequence import FastaSequence
		import matplotlib.pyplot as plt

		ascisse = []
		ascisse2 = []
		ascisse3 = []
		ascisse4 = []
		ascisse5 = []
		ordinate = []
		ordinate2 = []
		ordinate3 = []
		ordinate4 = []
		ordinate5 = []
		from main_estimate_jaccard import estimate_jaccard_difference_split

		# CLI inputs
		files = input_files()
		kmer_size = input("Tell me the kmer size (21) > ")
		kmer_size = 21 if kmer_size == "" else int(kmer_size)
		print()
		factorization = input("Tell me the factorization (cfl_icfl_comb)")
		if factorization == "":
			factorization = "cfl_icfl_comb"
		kfinger_size = input("tell me the kfinger sizes(3)")
		kfinger_size = 4 if kfinger_size == "" else int(kfinger_size)
		max = 0
		import copy

		y = 0

		files = sorted(files, key=lambda x: int(x.split(".")[0].split("/")[1]))

		for i in range(len(files) - 1):
			file1 = files[i]
			file2 = files[i + 1]
			print("{}%".format(i))
			# print("Comparing", file1, "-", file2)
			seq1 = FastaSequence(file1)
			seq2 = FastaSequence(file2)

			calc, estim = estimate_jaccard_difference_split(seq1, seq2, kmer_size, factorization, 3)
			# print(calc, estim)
			sketchsize = 1000

			import subprocess

			out = subprocess.Popen(["mash", "dist", file1, file2, "-s", str(sketchsize), "-k", str(kmer_size)],
			                       stdout=subprocess.PIPE,
			                       stderr=subprocess.STDOUT)
			stdout, stderr = out.communicate()
			# print(stdout)
			# output = str(stdout).split('t')[-1].split('\\')[0]
			# output = output.split("/")
			# mashresult1 = int(output[0]) / int(output[1])
			# print("mash kmer {}".format(mashresult1))
			output = str(stdout).split("/" + str(sketchsize))[0].split("t")[-1]
			mashresult1 = int(output) / sketchsize

			file1pre = file1.replace("/", "_preprocessed/preprocessed_")
			file2pre = file2.replace("/", "_preprocessed/preprocessed_")

			out = subprocess.Popen(
				["mash", "dist", file1pre, file2pre, "-s", str(sketchsize), "-k", str(kfinger_size), "-z",
				 "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZàèìòù%^-+:_çé!?£$§°*"],
				stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			stdout, stderr = out.communicate()
			# print(stdout)
			# print(stderr)
			output = str(stdout).split("/" + str(sketchsize))[0].split("t")[-1]
			try:
				mashresult2 = int(output) / sketchsize
			except ValueError:
				continue

			out = subprocess.Popen(
				["mash", "dist", file1pre, file2pre, "-s", str(sketchsize), "-k", str(kfinger_size + 1), "-z",
				 "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZàèìòù%^-+:_çé!?£$§°*"],
				stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			stdout, stderr = out.communicate()
			# print(stdout)
			# print(stderr)
			output = str(stdout).split('t')[-1].split('\\')[0]
			output = output.split("/")
			try:
				mashresult3 = int(output[0]) / int(output[1])
			except ValueError:
				continue

			out = subprocess.Popen(
				["mash", "dist", file1pre, file2pre, "-s", str(sketchsize), "-k", str(kfinger_size + 2), "-z",
				 "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZàèìòù%^-+:_çé!?£$§°*"],
				stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			stdout, stderr = out.communicate()
			# print(stdout)
			# print(stderr)
			output = str(stdout).split('t')[-1].split('\\')[0]
			output = output.split("/")
			try:
				mashresult4 = int(output[0]) / int(output[1])
			except ValueError:
				continue
			# print("mash kfing {}".format(mashresult2))

			# ordinate2.append(abs(mashresult1-mashresult2))
			# ascisse2.append(mashresult1)
			ordinate2.append(1 - abs(mashresult1 - calc))
			ascisse2.append(calc)

			ordinate3.append(1 - abs(mashresult2 - calc))
			ascisse3.append(calc)

			ordinate.append(1 - abs(estim - calc))
			ascisse.append(calc)

			ordinate4.append(1 - abs(mashresult3 - calc))
			ascisse4.append(calc)

			ordinate5.append(1 - abs(mashresult4 - calc))
			ascisse5.append(calc)

			print("---------------------------------")
			print()
			if abs(calc - estim) > 0.25:
				y += 1

		plt.figure(num=None, figsize=(15, 10), dpi=80, facecolor="w", edgecolor="k")
		# plt.plot(ascisse,	ordinate)
		plt.plot(ascisse2, ordinate2, "b.")
		plt.plot(ascisse, ordinate, "r+")

		plt.plot(ascisse3, ordinate3, "y*")
		plt.plot(ascisse4, ordinate4, "g*")
		plt.plot(ascisse5, ordinate5, "r*")
		# plt.axis([0, 1, 0, 1])
		plt.xlabel("Similarità stringhe")
		plt.ylabel("Precisione")
		plt.legend(["Mash " + str(kmer_size) + "-mer",
		            "Jaccard " + str(kfinger_size) + "-finger ",
		            "Mash " + str(kfinger_size) + "-finger",
		            "Mash " + str(kfinger_size + 1) + "-finger",
		            "Mash " + str(kfinger_size + 2) + "-finger"])
		plt.title("Figura A")
		plt.grid()
		plt.show()
		print("precisione medio con mash {}-mer = {}".format(str(kmer_size), sum(ordinate2) / len(ordinate2)))
		print("precisione medio con mash {}-finger = {}".format(str(kfinger_size), sum(ordinate3) / len(ordinate3)))
		print("precisione medio con mash {}-finger = {}".format(str(kfinger_size + 1), sum(ordinate4) / len(ordinate4)))
		print("precisione medio con mash {}-finger = {}".format(str(kfinger_size + 2), sum(ordinate5) / len(ordinate5)))
		print("precisione jaccard con {}-finger = {}".format(str(kfinger_size), sum(ordinate) / len(ordinate)))

	elif action == 9:
		def generate_parts(dir, maxlen):
			files = os.listdir(dir)
			files = [dir + "/" + file for file in files]

			cur_source_index = 0
			cur_source_data = None

			cur_dest_index = 0
			cur_dest_data = ""

			while cur_source_index < len(files):
				if cur_source_data is None:
					seq = FastaSequence(files[cur_source_index])
					cur_source_data = seq.get_data()
				data_to_consume_len = len(cur_source_data)
				if len(cur_dest_data) + data_to_consume_len <= maxlen:
					cur_dest_data += cur_source_data
					cur_source_data = None
					cur_source_index += 1
				else:
					space_to_fill = maxlen - len(cur_dest_data)
					cur_dest_data += cur_source_data[:space_to_fill]
					cur_source_data = cur_source_data[space_to_fill:]
					cur_dest_index += 1
					yield cur_dest_data
					cur_dest_data = ""
			if len(cur_dest_data) > 0:
				yield cur_dest_data


		dir = "prova"
		maxlen = 250
		for part in generate_parts(dir, maxlen):
			print(part)

	elif action == 10:

		alfabet = input("Tell me the alfabet (ACGT)> ")
		alfabet = list("ACGT") if alfabet == "" else list(alfabet)

		datasets_count = int(input("Tell me how much datasets do you want> "))
		datasets_folder = input("Tell me the dataset directory> ")

		print("Tell me the dataset size")
		print("  150-200:  short read")
		print("  5k-10k:   medium read")
		print("  50k-100k: long read")
		dataset_size = input("> ")
		dataset_size = int(dataset_size[:-1]) * 1000 if dataset_size.endswith("k") else int(dataset_size)

		import random
		import os

		if not os.path.isdir(datasets_folder):
			os.mkdir(datasets_folder)
		maxcharsonline = 70
		for dataset in range(1, datasets_count + 1):
			f = open(datasets_folder + "/" + str(dataset) + ".fasta", "w")
			f.write(">")
			i = 0
			while i < dataset_size:
				if i % maxcharsonline == 0:
					f.write("\n")
				f.write(random.choice(alfabet))
				i += 1
			f.write("\n")
			f.close()

	elif action == 11:

		from format import write_into_fasta
		from sequences.Sequence import FastaSequence

		option = input("[S] scrivi\n[K] Compara con k-mers\n[F] Compara con k-fingers\n> ")
		if option == "S":
			seq = FastaSequence("in/13.fasta")
			offset = 30
			width = 150
			for i in range(11):
				read = seq.get_data()[i * (width - offset): i * (width - offset) + width]
				write_into_fasta("dsa/" + str(i + 1) + ".fasta", read)

		elif option == "K":

			from suites.jaccard_kmers_suite import calc_similarity_with_kmers_from_fasta_pair, stdout_output_function

			kmer_size = 21
			for left_index in range(1, 11):
				fastapath1 = "dsa/" + str(left_index) + ".fasta"
				fastapath2 = "dsa/" + str(left_index + 1) + ".fasta"
				calc_similarity_with_kmers_from_fasta_pair(fastapath1, fastapath2, kmer_size, stdout_output_function)

		elif option == "F":
			from main_estimate_jaccard import estimate_jaccard
			from sequences.Sequence import FastaSequence

			seq1 = FastaSequence("dsa/1.fasta")
			seq2 = FastaSequence("dsa/2.fasta")
			kmer_size = 21
			factorization = "cfl_icfl_comb"

			for use_super_fp in [False, True]:
				for kfinger_size in range(2, 8 + 1):
					print("Use super fingerprint", use_super_fp)
					print("Kfinger size", kfinger_size)
					calc = jaccard_on_kmers(seq1.get_data(), seq2.get_data(), kmer_size)
					estim = estimate_jaccard(seq1, seq2, factorization, kfinger_size, use_super_fp)
					print("calc", calc)
					print("estim", estim)
					print()

	else:
		print("Goodbye...️")
