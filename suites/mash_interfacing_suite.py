import os
from factorization import get_factors
from files_manager import input_files
from lyndon.utils import _complement
from sequences.Sequence import FastaSequence
from traditional_technique import jaccard_on_kmers
import matplotlib.pyplot as plt
import subprocess
from main_estimate_jaccard import estimate_jaccard_difference_split


def mash_interfacing_suite_main(action):
	if action == 5:
		preprocess_directory()
	elif action == 7:
		grafico_mash_on_kmers_and_kfingers_on_preprocessed_dataset()
	elif action == 8:
		grafico_mash_and_jaccard_with_varying_kfingers()
	else:
		return False


def preprocess_directory():
	directory = input("Tell me the directory with the file to preprocess> ")

	factorization = "cfl_icfl_comb"
	# input("Tell me the algorithm for factoring\n   (cfl,icfl,cfl_icfl,cfl_comb,icfl_comb,cfl_icfl_comb)> ")

	use_super_fp = input("Tell me if you want to use the super-fingerprint [Y/n]> ")
	use_super_fp = (use_super_fp == "Y" or use_super_fp == "y" or use_super_fp == "")

	alf = "#0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZàèìòù%^-+:_çé!?£$§°*"


	files = os.listdir(directory)
	files.sort()
	for file in files:
		if file.endswith(".fasta"):
			print(file)
			seq = FastaSequence(directory + "/" + file)
			fingerprint = [len(f) for f in get_factors(factorization, seq.get_data())]
			if use_super_fp:
				fingerprint += [0]
				fingerprint += [len(f) for f in get_factors(factorization, _complement(seq.get_data()))]
			if max(fingerprint) >= len(alf):
				print('   la fingerprint non può essere codificata mediante la ristretta porzione di ASCII predisposta')
				continue
			if not os.path.isdir(directory + "_preprocessed"):
				os.mkdir(directory + "_preprocessed")
			savefile = open(directory + "_preprocessed/preprocessed_" + file, "w")
			savefile.write(">\'" + seq.get_title() + "\' preprocessed encoding Lyndon ")
			if use_super_fp:
				savefile.write("super-")
			savefile.write("fingerprint with ASCII code\n")
			for finger in fingerprint:
				savefile.write(alf[finger])
			savefile.write('\n')
			savefile.close()


def grafico_mash_on_kmers_and_kfingers_on_preprocessed_dataset():
	files = input_files()

	kmer_size = input("Tell me the kmer size (21) > ")
	kmer_size = 21 if kmer_size == "" else int(kmer_size)

	factorization = input("Tell me the factorization (cfl_icfl_comb)")
	if factorization == "":
		factorization = "cfl_icfl_comb"

	kfinger_size = input("Tell me the kfinger sizes(3)")
	kfinger_size = 3 if kfinger_size == "" else int(kfinger_size)

	# total_max = 0
	# for file in files:
	# 	seq = FastaSequence(file)
	# 	read = seq.get_data()
	# 	factors = [len(f) for f in get_factors(factorization, read)]
	# 	factors2 = [len(f) for f in get_factors(factorization, _complement(read))]
	# 	current_max = max(max(factors), max(factors2))
	# 	if current_max > total_max:
	# 		total_max = current_max
	# print("Il fattore più grande ha lunghezza {}".format(total_max))
	# TODO: o metti nel preprocessed o no (Antonio ha fatto in modo che i file che preprocessiamo erano gia con lunghezze giuste, quindi non dava mai problemi)
	# Generalmente, il preprocessing non va proprio a buon fine, se i fattori sono troppo grandi...

	files = sorted(files, key=lambda x: int(x.split(".")[0].split("/")[1]))

	x_mash_on_kmers = []
	y_mash_on_kmers = []

	# x_jaccard_on_kfingers = []
	# y_jaccard_on_kfingers = []

	x_mash_on_kfingers = []
	y_mash_on_kfingers = []

	for i in range(len(files) - 1):
		file1 = files[i]
		file2 = files[i + 1]

		print("Comparing", file1, "-", file2)
		seq1 = FastaSequence(file1)
		seq2 = FastaSequence(file2)

		calc = jaccard_on_kmers(seq1.get_data(), seq2.get_data(), kmer_size)

		# Mash on k-mers
		out = subprocess.Popen(["mash", "dist", file1, file2, "-s", "1000"], stdout=subprocess.PIPE,
		                       stderr=subprocess.STDOUT)
		stdout, stderr = out.communicate()
		output = str(stdout).split('t')[-1].split('\\')[0]
		output = output.split("/")
		mashresult1 = int(output[0]) / int(output[1])
		x_mash_on_kmers.append(calc)
		y_mash_on_kmers.append(1 - abs(calc - mashresult1))

		# Mash on k-fingers
		file1pre = file1.replace("/", "_preprocessed/preprocessed_")
		file2pre = file2.replace("/", "_preprocessed/preprocessed_")
		out = subprocess.Popen(["mash", "dist", file1pre, file2pre, "-s", "1000", "-k", "3", "-z",
		                        "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZàèìòù%^-+:_çé!?£$§°*"],
		                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		stdout, stderr = out.communicate()
		output = str(stdout).split('t')[-1].split('\\')[0]
		output = output.split("/")
		mashresult2 = int(output[0]) / int(output[1])
		x_mash_on_kfingers.append(calc)
		y_mash_on_kfingers.append(1 - abs(calc - mashresult2))

	plt.title("Factorization: " + factorization + "; kfinger size: " + str(kfinger_size))
	plt.axis([0, 1, 0, 1])
	plt.xlabel("Similarità stringhe")
	plt.ylabel("Similarità stimata")
	plt.plot(x_mash_on_kmers, y_mash_on_kmers, "b.")
	plt.plot(x_mash_on_kfingers, y_mash_on_kfingers, "r.")
	plt.legend(["Mash k-mer", "Mash k-finger"])
	plt.grid()
	plt.show()


def grafico_mash_and_jaccard_with_varying_kfingers():
	x_jaccard_on_4fingers = []
	x_mash_on_kmers = []
	x_mash_on_4fingers = []
	x_mash_on_5fingers = []
	x_mash_on_6fingers = []

	y_jaccard_on_4fingers = []
	y_mash_on_kmers = []
	y_mash_on_4fingers = []
	y_mash_on_5fingers = []
	y_mash_on_6fingers = []

	files = input_files()

	kmer_size = input("Tell me the kmer size (21) > ")
	kmer_size = 21 if kmer_size == "" else int(kmer_size)

	factorization = input("Tell me the factorization (cfl_icfl_comb)")  # TODO: > si o no?
	if factorization == "":
		factorization = "cfl_icfl_comb"

	kfinger_size = input("Tell me the kfinger sizes (3)")
	kfinger_size = 4 if kfinger_size == "" else int(kfinger_size)

	# TODO: voleva sapere il total_max pure quì

	sketchsize = 1000
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

		out = subprocess.Popen(["mash", "dist", file1, file2, "-s", str(sketchsize), "-k", str(kmer_size)],
		                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
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

		# x_mash_on_kmers.append(mashresult1)
		# y_mash_on_kmers.append(abs(mashresult1-mashresult2))

		x_mash_on_kmers.append(calc)
		y_mash_on_kmers.append(1 - abs(mashresult1 - calc))

		x_mash_on_4fingers.append(calc)
		y_mash_on_4fingers.append(1 - abs(mashresult2 - calc))

		x_jaccard_on_4fingers.append(calc)
		y_jaccard_on_4fingers.append(1 - abs(estim - calc))

		x_mash_on_5fingers.append(calc)
		y_mash_on_5fingers.append(1 - abs(mashresult3 - calc))

		x_mash_on_6fingers.append(calc)
		y_mash_on_6fingers.append(1 - abs(mashresult4 - calc))

	plt.figure(num=None, figsize=(15, 10), dpi=80, facecolor="w", edgecolor="k")
	plt.title("Factorization: " + factorization + "; kfinger size: " + str(kfinger_size) + "; sketch size: " + str(
		sketchsize))
	# plt.axis([0, 1, 0, 1])
	plt.xlabel("Similarità stringhe")
	plt.ylabel("Precisione")
	plt.plot(x_mash_on_kmers, y_mash_on_kmers, "b.")
	plt.plot(x_jaccard_on_4fingers, y_jaccard_on_4fingers, "r+")
	plt.plot(x_mash_on_4fingers, y_mash_on_4fingers, "y*")
	plt.plot(x_mash_on_5fingers, y_mash_on_5fingers, "g*")
	plt.plot(x_mash_on_6fingers, y_mash_on_6fingers, "r*")
	plt.legend(["Mash " + str(kmer_size) + "-mer",
	            "Jaccard " + str(kfinger_size) + "-finger ",
	            "Mash " + str(kfinger_size) + "-finger",
	            "Mash " + str(kfinger_size + 1) + "-finger",
	            "Mash " + str(kfinger_size + 2) + "-finger"])
	plt.grid()
	plt.show()

	print("\nPrecisione medio con mash %-2d-mers:    %7.5f" % (kmer_size, sum(y_mash_on_kmers) / len(y_mash_on_kmers)))
	print("Precisione jaccard con %-2d-fingers:    %7.5f" % (
		kfinger_size, sum(y_jaccard_on_4fingers) / len(y_jaccard_on_4fingers)))
	print("Precisione medio con mash %-2d-fingers: %7.5f" % (
		kfinger_size, sum(y_mash_on_4fingers) / len(y_mash_on_4fingers)))
	print("Precisione medio con mash %-2d-fingers: %7.5f" % (
		kfinger_size + 1, sum(y_mash_on_5fingers) / len(y_mash_on_5fingers)))
	print("Precisione medio con mash %-2d-fingers: %7.5f" % (
		kfinger_size + 2, sum(y_mash_on_6fingers) / len(y_mash_on_6fingers)))