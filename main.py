#!/bin/env python3

import main_calculate_jaccard
import main_export_fq_in_fasta_files
import main_mash_interface
import main_process_files_in_directory

if __name__ == "__main__":
	print("1: Export FASTA from FQ")
	print("2: Calculate Jaccard similarity (exploring various kmer size)")
	print("3: Calculate Jaccard similarity (fixed kmer size)")
	print("4: Estimate Jaccard similarity (estimate on pairs of some files)")
	print("5: Preprocess all the files in a directory")
	print("6: Encode in ASCII the fingers of two FASTA files")
	print("7: WE (pure a me plz)")
	print("8: Antonio dammi un nome...")
	print("9: Crea dataset di long long reads")
	print("_: Quit")
	action = int(input("\nTell me what to do > "))
	if action == 1:
		main_export_fq_in_fasta_files.main()
	elif action == 2:
		main_calculate_jaccard.explore_kmer_size()
	elif action == 3:
		main_calculate_jaccard.fixed_kmer_size()
	elif action == 4:

		from main_estimate_jaccard import jaccard_thanks_factorizations, get_csv_exporter
		import os
		from collections import defaultdict
		from itertools import combinations
		from Sequence import FastaSequence

		# CLI inputs
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

		kmer_size = input("Tell me the kmer size (21) > ")
		kmer_size = 21 if kmer_size == "" else int(kmer_size)
		tolerance = input("Tell me the tolerance (0.15) > ")
		tolerance = 0.15 if tolerance == "" else float(tolerance)
		print()

		# Creation of a work directory
		dirname = []
		for path in files:
			full_filename = path.split("/")[-1]
			filename_parts_no_extension = path.split("/")[-1].split(".")[:-1]
			dirname.append(".".join(filename_parts_no_extension))
		dirname = "-".join(dirname) + " kmer " + str(kmer_size) + "-tolerance " + str(int(tolerance * 100))
		attempt = 1
		if os.path.isdir(dirname):
			attempt = 2
			dirname = dirname + " (2)"
		while os.path.isdir(dirname):
			attempt += 1
			dirname = dirname[:dirname.rfind("(")] + "(" + str(attempt) + ")"
		os.mkdir(dirname)

		# Initializing the CSV export
		csv_exporter = get_csv_exporter(dirname)

		summary = open(dirname + "/summary.txt", "w")

		sopravvissuti = defaultdict(int)
		for file1, file2 in combinations(files, 2):
			print("Comparing", file1, "-", file2)

			seq1 = FastaSequence(file1)
			seq2 = FastaSequence(file2)

			result = jaccard_thanks_factorizations(seq1, seq2, kmer_size, tolerance, csv_exporter)

			for factorization, configs in result["factorizations"].items():
				for split, window_size in configs:
					sopravvissuti[(split, window_size, factorization)] += 1

			summary.write(seq1.get_name() + " - " + seq2.get_name() + "\n")
			summary.write(seq1.get_title() + "\n")
			summary.write(seq2.get_title() + "\n\n")
			summary.write("kmer size ->    " + str(kmer_size) + "\n")
			summary.write("jaccard   ->    %5.2f" % (result["jaccard"] * 100) + "%\n")
			summary.write("tolerance ->    %5.2f" % (tolerance * 100) + "%\n\n")
			for factorization, configs in result["factorizations"].items():
				summary.write("  %13s  %3d\n" % (factorization, len(configs)))
			summary.write("\n----------------------------------------------------------------------\n\n")

		sopravvissuti = sorted(sopravvissuti.items(), key=lambda p: p[1], reverse=True)
		for (split, window_size, factorization), count in sopravvissuti:
			summary.write("%3d-%1d %15s   -> %3d\n" % (split, window_size, factorization, count))
		summary.close()

	elif action == 5:
		main_process_files_in_directory.main()
	elif action == 6:
		main_mash_interface.main()
	elif action == 7:
		from main_estimate_jaccard import jaccard_thanks_factorizations, get_csv_exporter
		import os
		from collections import defaultdict
		from itertools import combinations
		from Sequence import FastaSequence
		import matplotlib.pyplot as plt

		ascisse = []
		ascisse2 = []
		ascisse3 = []
		ordinate = []
		ordinate2 = []
		ordinate3 = []
		from main_estimate_jaccard import estimate_jaccard_difference

		# CLI inputs
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

		kmer_size = input("Tell me the kmer size (21) > ")
		kmer_size = 21 if kmer_size == "" else int(kmer_size)
		print()
		factorization = input("Tell me the factorization (cfl_icfl_comb)")
		if factorization == "":
			factorization = "cfl_icfl_comb"
		kfinger_size = input("tell me the kfinger sizes(3)")
		kfinger_size = 3 if kfinger_size == "" else int(kfinger_size)
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

			calc, estim = estimate_jaccard_difference(seq1, seq2, kmer_size, factorization, kfinger_size)
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
		from factorization import get_fingers_after_split
		from main_estimate_jaccard import jaccard_thanks_factorizations, get_csv_exporter
		import os
		from collections import defaultdict
		from itertools import combinations
		from Sequence import FastaSequence
		import matplotlib.pyplot as plt

		ascisse = []
		ascisse2 = []
		ascisse3 = []
		ordinate = []
		ordinate2 = []
		ordinate3 = []
		from main_estimate_jaccard import estimate_jaccard_difference_split

		# CLI inputs
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

		kmer_size = input("Tell me the kmer size (21) > ")
		kmer_size = 21 if kmer_size == "" else int(kmer_size)
		print()
		factorization = input("Tell me the factorization (cfl_icfl_comb)")
		if factorization == "":
			factorization = "cfl_icfl_comb"
		kfinger_size = input("tell me the kfinger sizes(3)")
		kfinger_size = 3 if kfinger_size == "" else int(kfinger_size)
		max = 0
		import copy

		y = 0

		files = sorted(files, key=lambda x: int(x.split(".")[0].split("/")[1]))

		for i in range(len(files) - 1):
			file1 = files[i]
			file2 = files[i + 1]

			print("Comparing", file1, "-", file2)
			seq1 = FastaSequence(file1)
			seq2 = FastaSequence(file2)

			calc, estim = estimate_jaccard_difference_split(seq1, seq2, kmer_size, factorization, kfinger_size)
			print(calc, estim)

			import subprocess

			out = subprocess.Popen(["mash", "dist", file1, file2, "-s", "10"], stdout=subprocess.PIPE,
			                       stderr=subprocess.STDOUT)
			stdout, stderr = out.communicate()
			print(stdout)
			output = str(stdout).split('t')[-1].split('\\')[0]
			output = output.split("/")
			mashresult1 = int(output[0]) / int(output[1])
			print("mash kmer {}".format(mashresult1))

			file1pre = file1.replace("/", "_preprocessed/preprocessed_")
			file2pre = file2.replace("/", "_preprocessed/preprocessed_")

			out = subprocess.Popen(["mash", "dist", file1pre, file2pre, "-s", "10", "-k", "4", "-z",
			                        "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZàèìòù%^-+:_çé!?£$§°*"],
			                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			stdout, stderr = out.communicate()

			output = str(stdout).split('t')[-1].split('\\')[0]
			output = output.split("/")
			mashresult2 = int(output[0]) / int(output[1])
			print("mash kfing {}".format(mashresult2))

			# ordinate2.append(abs(mashresult1-mashresult2))
			# ascisse2.append(mashresult1)
			ordinate2.append(1 - abs(mashresult1 - calc))
			ascisse2.append(calc)

			ordinate3.append(1 - abs(mashresult2 - calc))
			ascisse3.append(calc)

			ordinate.append(1 - abs(estim - calc))
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

	elif action == 9:

		import m_create_longlongreads

		dir = "prova"
		maxlen = 250
		for part in m_create_longlongreads.main(dir, maxlen):
			print(part)

	else:
		print("Goodbye...️")
