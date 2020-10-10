#!/bin/env python3

from files_manager import input_files_in_directory, input_factorization, input_files
from suites.jaccard_kfingers_suite import great_estimation
from suites.jaccard_kmers_suite import progressing_jaccard_on_kmers, stdout_output
from suites.mash_interfacing_suite import preprocess_directory, \
	show_plot_mash_on_kmers_and_kfingers_based_on_preprocessed_dataset, \
	show_plot_mash_jaccard_on_varying_kfingers_preprocessed_dataset
from traditional_technique import monitor_jaccard_on_kmers_with_sequences, stdout_jaccard_on_kmers_output_function
from sequences.Sequence import FastaSequence
from itertools import combinations

if __name__ == '__main__':
	print('JACCARD ON K-MERS OR K-FINGERS')
	print('1: Calculate Jaccard similarity (exploring various k-mer size)')
	print('2: Calculate Jaccard similarity (fixed k-mer size)')
	print('3: Estimate Jaccard similarity (estimate on pairs of some files)')
	print('')
	print('MASH INTERFACING')
	print('4: Preprocess all the files in a directory (only short reads)')
	print('5: Graph with MASH on k-mers and k-fingers (on preprocessed dataset) for short reads')
	print('6: Graph with MASH (and Jaccard) on k-mers and varying k-fingers (on preprocessed dataset) for long reads')
	print('')
	print('_: Quit')
	print('')
	action = int(input('\nTell me what to do> '))

	if action == 1:
		seq1 = FastaSequence(input('Tell me the first file> '))
		seq2 = FastaSequence(input('Tell me the second file> '))
		print('Comparing:', seq1.get_name(), 'with', seq2.get_name(), '\n')
		progressing_jaccard_on_kmers(seq1.get_data(), seq2.get_data(), 10, stdout_output)

	elif action == 2:
		files = input_files()
		kmer_size = input('Tell me the k-mer size (usually 21)> ')
		kmer_size = 21 if kmer_size == '' else int(kmer_size)
		sequences = [FastaSequence(file) for file in files]
		for seq1, seq2 in combinations(sequences, 2):
			monitor_jaccard_on_kmers_with_sequences(seq1, seq2, kmer_size, stdout_jaccard_on_kmers_output_function)

	elif action == 3:
		files = input_files()
		kmer_size = input('Tell me the k-mer size (21)> ')
		kmer_size = 21 if kmer_size == '' else int(kmer_size)
		tolerance = input('Tell me the tolerance (0.15)> ')
		tolerance = 0.15 if tolerance == '' else float(tolerance)
		great_estimation(files, kmer_size, tolerance)

	elif action == 4:
		files = input_files_in_directory('Tell me the directory with the file to preprocess')
		factorization = input_factorization()
		use_super_fp = input('Tell me if you want to use the super-fingerprint [Y/n]> ')
		use_super_fp = use_super_fp == 'Y' or use_super_fp == 'y' or use_super_fp == ''
		preprocess_directory(files, factorization, use_super_fp)

	elif action == 5:
		files = input_files_in_directory('Tell me the dataset directory (that corresponds to the preprocessed dataset)')
		kmer_size = input('Tell me the kmer size (21)> ')
		kmer_size = 21 if kmer_size == '' else int(kmer_size)
		show_plot_mash_on_kmers_and_kfingers_based_on_preprocessed_dataset(files, kmer_size)

	elif action == 6:
		files = input_files_in_directory('Tell me the dataset directory (that corresponds to the preprocessed dataset)')
		kmer_size = input('Tell me the kmer size (21)> ')
		kmer_size = 21 if kmer_size == '' else int(kmer_size)
		factorization = input_factorization()
		kfinger_size = input('Tell me the kfinger sizes to use on MASH (4)> ')
		kfinger_size = 4 if kfinger_size == '' else int(kfinger_size)
		show_plot_mash_jaccard_on_varying_kfingers_preprocessed_dataset(files, kmer_size, factorization, kfinger_size)

	else:
		print('Goodbye...️')
