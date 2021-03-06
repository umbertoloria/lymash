import os
import subprocess

import matplotlib.pyplot as plt

from new_combined_technique import create_fingerprint, jaccard_on_kfingers
from sequences import FastaSequence
from traditional_technique import jaccard_on_kmers


def preprocess_directory(files: list, factorization: str, use_super_fp: bool, split: int):
	"""Creates a preprocessed directory containing preprocessed files related to the specified files. A preprocessed
	file is a (improper, since the alphabet is not the standard one) FASTA file containig an encoded fingerprint of the
	related normal file. The model of the used fingerprint is defined by the 'factorization', 'use_super_fp' and 'split'
	parameters."""

	if len(files) == 0:
		print('The directory can\'t be empty')
		return
	preprocessed_dataset_directory = os.path.dirname(files[0]) + '_preprocessed'
	if not os.path.isdir(preprocessed_dataset_directory):
		os.mkdir(preprocessed_dataset_directory)

	for file in files:
		print(file)
		seq = FastaSequence(file)
		prepr = preprocess_string(seq.get_data(), factorization, use_super_fp, split)
		if prepr is not None:
			savefile = open('_preprocessed/preprocessed_'.join(os.path.split(file)), 'w')
			savefile.write('>\'' + seq.get_title() + '\' preprocessed encoding Lyndon ')
			if use_super_fp:
				savefile.write('super-')
			savefile.write('fingerprint (' + factorization + ') with ASCII code\n')
			savefile.write(prepr)
			savefile.write('\n')
			savefile.close()
		else:
			print('   La fingerprint non può essere codificata mediante la ristretta porzione di ASCII predisposta')


def preprocess_string(text: str, factorization: str, use_super_fp: bool, split: int):
	alf = '#0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZàèìòù%^-+:_çé!?£$§°*'
	fingerprint = create_fingerprint(text, factorization, use_super_fp, split)
	if max(fingerprint) < len(alf):
		return ''.join([alf[f] for f in fingerprint])
	else:
		return None


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
	output = stdout.split('t')[-1].split('\\')[0].split('/')
	mashresult = int(output[0]) / int(output[1])
	return mashresult


def show_plot_mash_on_kmers_and_kfingers_based_on_preprocessed_dataset(files: list, kmer_size: int):
	"""This function creates a plot containing comparative results using MASH on normal files (dividing them in k-mers)
	and preprocessed files (dividing them in k-mers of an encoded fingerprint, actually k-fingers). The 'files'
	parameter contains filepaths of the normal files (of which already exists a preprocessed counterpart). The length of
	the (actual) k-mers is defined by the 'kmer_size' parameter. This function should be used on an already preprocessed
	short reads dataset only."""

	x_mash_on_kmers = []
	y_mash_on_kmers = []

	x_mash_on_kfingers = []
	y_mash_on_kfingers = []

	alf = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZàèìòù%^-+:_çé!?£$§°*'

	for i in range(len(files) - 1):
		file1 = files[i]
		file2 = files[i + 1]

		print('Comparing', file1, '<->', file2)
		seq1 = FastaSequence(file1)
		seq2 = FastaSequence(file2)

		calc = jaccard_on_kmers(seq1.get_data(), seq2.get_data(), kmer_size)

		# Mash on k-mers
		x_mash_on_kmers.append(calc)
		mashresult1 = use_mash(file1, file2, 1000)
		y_mash_on_kmers.append(1 - abs(calc - mashresult1))

		# Mash on k-fingers
		file1pre = '_preprocessed/preprocessed_'.join(os.path.split(file1))
		file2pre = '_preprocessed/preprocessed_'.join(os.path.split(file2))
		x_mash_on_kfingers.append(calc)
		mashresult2 = use_mash(file1pre, file2pre, 1000, 3, alf)
		y_mash_on_kfingers.append(1 - abs(calc - mashresult2))

	plt.title('Comparison with k-mer size: ' + str(kmer_size))
	plt.axis([0, 1, 0, 1.05])
	plt.xlabel('Similarità stringhe')
	plt.ylabel('Similarità stimata')
	plt.plot(x_mash_on_kmers, y_mash_on_kmers, 'b.')
	plt.plot(x_mash_on_kfingers, y_mash_on_kfingers, 'r.')
	plt.legend(['Mash k-mer', 'Mash k-finger'])
	plt.grid()
	plt.show()


def show_plot_mash_jaccard_varying_kfingers(files: list, kmer_size: int, factorization: str, kfinger_size: int,
                                            use_super_fp: bool, split: int):
	"""Very similar to the 'show_plot_mash_on_kmers_and_kfingers_based_on_preprocessed_dataset' function, but unlike
	that, this should be used on an already preprocessed long reads dataset. Multiple k-fingers length will be used
	(precisely, kfinger_size, kfinger_size+1 and kfinger_size+2). In addition to the MASH single k-mers and multiple
	k-fingers usages, it will be also calculated the Jaccard similarity on k-fingers (on size 'kfinger_size' only)."""

	alf = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZàèìòù%^-+:_çé!?£$§°*'
	sketchsize = 1000

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
		print('Comparing', file1, '<->', file2)
		seq1 = FastaSequence(file1)
		seq2 = FastaSequence(file2)

		calc = jaccard_on_kmers(seq1.get_data(), seq2.get_data(), kmer_size)
		estim = jaccard_on_kfingers(seq1.get_data(), seq2.get_data(), factorization, 3, use_super_fp, split)
		x_jaccard_on_3fingers.append(calc)
		y_jaccard_on_3fingers.append(1 - abs(estim - calc))

		mashresult1 = use_mash(file1, file2, sketchsize, kmer_size)
		x_mash_on_kmers.append(calc)
		y_mash_on_kmers.append(1 - abs(mashresult1 - calc))

		file1pre = '_preprocessed/preprocessed_'.join(os.path.split(file1))
		file2pre = '_preprocessed/preprocessed_'.join(os.path.split(file2))

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

	plt.figure(num=None, figsize=(15, 10), dpi=80, facecolor='w', edgecolor='k')
	plt.title('Factorization: ' + factorization + '; kfinger size: ' + str(kfinger_size) + '; sketch size: ' + str(
		sketchsize))
	plt.xlabel('Similarità stringhe')
	plt.ylabel('Precisione')
	plt.plot(x_mash_on_kmers, y_mash_on_kmers, 'b.')
	plt.plot(x_jaccard_on_3fingers, y_jaccard_on_3fingers, 'r+')
	plt.plot(x_mash_on_4fingers, y_mash_on_4fingers, 'y*')
	plt.plot(x_mash_on_5fingers, y_mash_on_5fingers, 'g*')
	plt.plot(x_mash_on_6fingers, y_mash_on_6fingers, 'r*')
	legend_list = ['Mash ' + str(kmer_size) + '-mer', 'Jaccard 3-finger', 'Mash ' + str(kfinger_size) + '-finger']
	legend_list += ['Mash ' + str(kfinger_size + 1) + '-finger', 'Mash ' + str(kfinger_size + 2) + '-finger']
	plt.legend(legend_list)
	plt.grid()
	plt.show()

	print('\nPrecisione medio con mash %d-mers: %7.5f' % (kmer_size, sum(y_mash_on_kmers) / len(y_mash_on_kmers)))
	print('Precisione jaccard con %d-fingers: %7.5f' % (
		kfinger_size, sum(y_jaccard_on_3fingers) / len(y_jaccard_on_3fingers)))
	print('Precisione medio con mash %d-fingers: %7.5f' % (
		kfinger_size, sum(y_mash_on_4fingers) / len(y_mash_on_4fingers)))
	print('Precisione medio con mash %d-fingers: %7.5f' % (
		kfinger_size + 1, sum(y_mash_on_5fingers) / len(y_mash_on_5fingers)))
	print('Precisione medio con mash %d-fingers: %7.5f' % (
		kfinger_size + 2, sum(y_mash_on_6fingers) / len(y_mash_on_6fingers)))
