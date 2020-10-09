import os
from factorization import get_factors
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


def use_mash(filepath1: str, filepath2: str, sketch_size: int = 0, window_size: int = 0, alphabet: str = ''):
	args_of_call = ['mash', 'dist', filepath1, filepath2]
	if sketch_size > 0:
		args_of_call += ['-s', str(sketch_size)]
	if window_size > 0:
		args_of_call += ['-k', str(window_size)]
	if alphabet != '':
		args_of_call += ['-z', alphabet]
	out = subprocess.Popen(args_of_call, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	stdout, stderr = out.communicate()
	stdout = str(stdout)
	if 'ERROR' in stdout:
		return -1
	if 'WARNING' in stdout:
		stdout = stdout[0:stdout.rfind('WARNING') - 1]
	output = stdout.split('t')[-1].split('\\')[0]
	output = output.split("/")
	mashresult = int(output[0]) / int(output[1])
	return mashresult


def grafico_mash_on_kmers_and_kfingers_on_preprocessed_dataset():
	directory = input("Tell me the directory with normal files (corresponding preprocessed files will be figured)> ")
	files = [directory + '/' + filename for filename in os.listdir(directory)]
	files.sort()

	kmer_size = input("Tell me the kmer size (21)> ")
	kmer_size = 21 if kmer_size == "" else int(kmer_size)

	factorization = input("Tell me the factorization (cfl_icfl_comb)> ")
	if factorization == "":
		factorization = "cfl_icfl_comb"

	kfinger_size = input("Tell me the kfinger size (3)> ")
	kfinger_size = 3 if kfinger_size == "" else int(kfinger_size)

	x_mash_on_kmers = []
	y_mash_on_kmers = []

	x_mash_on_kfingers = []
	y_mash_on_kfingers = []

	alf = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZàèìòù%^-+:_çé!?£$§°*"

	for i in range(len(files) - 1):
		file1 = files[i]
		file2 = files[i + 1]

		print("Comparing", file1, "-", file2)
		seq1 = FastaSequence(file1)
		seq2 = FastaSequence(file2)

		calc = jaccard_on_kmers(seq1.get_data(), seq2.get_data(), kmer_size)

		# Mash on k-mers
		x_mash_on_kmers.append(calc)
		mashresult1 = use_mash(file1, file2, 1000)
		y_mash_on_kmers.append(1 - abs(calc - mashresult1))

		# Mash on k-fingers
		file1pre = file1[0:file1.rfind('/')] + '_preprocessed/preprocessed_' + file1[file1.rfind('/') + 1:]
		file2pre = file2[0:file2.rfind('/')] + '_preprocessed/preprocessed_' + file2[file2.rfind('/') + 1:]
		x_mash_on_kfingers.append(calc)
		mashresult2 = use_mash(file1pre, file2pre, 1000, 3, alf)
		y_mash_on_kfingers.append(1 - abs(calc - mashresult2))

	plt.title("Factorization: " + factorization + "; kfinger size: " + str(kfinger_size))
	plt.axis([0, 1, 0, 1.05])
	plt.xlabel("Similarità stringhe")
	plt.ylabel("Similarità stimata")
	plt.plot(x_mash_on_kmers, y_mash_on_kmers, "b.")
	plt.plot(x_mash_on_kfingers, y_mash_on_kfingers, "r.")
	plt.legend(["Mash k-mer", "Mash k-finger"])
	plt.grid()
	plt.show()


def grafico_mash_and_jaccard_with_varying_kfingers():
	directory = input("Tell me the directory with normal files (corresponding preprocessed files will be figured)> ")
	files = [directory + '/' + filename for filename in os.listdir(directory)]
	files.sort()

	kmer_size = input("Tell me the kmer size (21)> ")
	kmer_size = 21 if kmer_size == "" else int(kmer_size)

	factorization = input("Tell me the factorization (cfl_icfl_comb)> ")
	if factorization == "":
		factorization = "cfl_icfl_comb"

	kfinger_size = input("Tell me the kfinger sizes to use on MASH (4)> ")
	kfinger_size = 4 if kfinger_size == "" else int(kfinger_size)

	sketchsize = 1000

	alf = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZàèìòù%^-+:_çé!?£$§°*"

	x_mash_on_kmers = []
	y_mash_on_kmers = []

	x_jaccard_on_3fingers = []
	y_jaccard_on_3fingers = []

	x_mash_on_4fingers = []
	y_mash_on_4fingers = []

	x_mash_on_5fingers = []
	y_mash_on_5fingers = []

	x_mash_on_6fingers = []
	y_mash_on_6fingers = []

	for i in range(len(files) - 1):
		file1 = files[i]
		file2 = files[i + 1]
		print("Comparing", file1, "-", file2)
		seq1 = FastaSequence(file1)
		seq2 = FastaSequence(file2)

		calc, estim = estimate_jaccard_difference_split(seq1, seq2, kmer_size, factorization, 3)
		x_jaccard_on_3fingers.append(calc)
		y_jaccard_on_3fingers.append(1 - abs(estim - calc))

		mashresult1 = use_mash(file1, file2, sketchsize, kmer_size)
		x_mash_on_kmers.append(calc)
		y_mash_on_kmers.append(1 - abs(mashresult1 - calc))

		file1pre = file1.replace("/", "_preprocessed/preprocessed_")
		file2pre = file2.replace("/", "_preprocessed/preprocessed_")

		mashresult2 = use_mash(file1pre, file2pre, sketchsize, kfinger_size, alf)
		if mashresult2 == -1:
			continue
		x_mash_on_4fingers.append(calc)
		y_mash_on_4fingers.append(1 - abs(mashresult2 - calc))

		mashresult3 = use_mash(file1pre, file2pre, sketchsize, kfinger_size + 1, alf)
		if mashresult3 == -1:
			continue
		x_mash_on_5fingers.append(calc)
		y_mash_on_5fingers.append(1 - abs(mashresult3 - calc))

		mashresult4 = use_mash(file1pre, file2pre, sketchsize, kfinger_size + 2, alf)
		if mashresult4 == -1:
			continue
		x_mash_on_6fingers.append(calc)
		y_mash_on_6fingers.append(1 - abs(mashresult4 - calc))

	plt.figure(num=None, figsize=(15, 10), dpi=80, facecolor="w", edgecolor="k")
	plt.title("Factorization: " + factorization + "; kfinger size: " + str(kfinger_size) + "; sketch size: " + str(
		sketchsize))
	plt.xlabel("Similarità stringhe")
	plt.ylabel("Precisione")
	plt.plot(x_mash_on_kmers, y_mash_on_kmers, "b.")
	plt.plot(x_jaccard_on_3fingers, y_jaccard_on_3fingers, "r+")
	plt.plot(x_mash_on_4fingers, y_mash_on_4fingers, "y*")
	plt.plot(x_mash_on_5fingers, y_mash_on_5fingers, "g*")
	plt.plot(x_mash_on_6fingers, y_mash_on_6fingers, "r*")
	legend_list = ["Mash " + str(kmer_size) + "-mer", "Jaccard 3-finger", "Mash " + str(kfinger_size) + "-finger"]
	legend_list += ["Mash " + str(kfinger_size + 1) + "-finger", "Mash " + str(kfinger_size + 2) + "-finger"]
	plt.legend(legend_list)
	plt.grid()
	plt.show()

	print("\nPrecisione medio con mash %d-mers: %7.5f" % (kmer_size, sum(y_mash_on_kmers) / len(y_mash_on_kmers)))
	print("Precisione jaccard con %d-fingers: %7.5f" % (
		kfinger_size, sum(y_jaccard_on_3fingers) / len(y_jaccard_on_3fingers)))
	print("Precisione medio con mash %d-fingers: %7.5f" % (
		kfinger_size, sum(y_mash_on_4fingers) / len(y_mash_on_4fingers)))
	print("Precisione medio con mash %d-fingers: %7.5f" % (
		kfinger_size + 1, sum(y_mash_on_5fingers) / len(y_mash_on_5fingers)))
	print("Precisione medio con mash %d-fingers: %7.5f" % (
		kfinger_size + 2, sum(y_mash_on_6fingers) / len(y_mash_on_6fingers)))
