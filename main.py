#!/bin/env python3

from cli_input import *
from std_cli_inputs import *
from suites.dataset_creation import *
from suites.jaccard_kfingers_suite import *
from suites.jaccard_kmers_suite import *
from suites.mash_interfacing_suite import *
from traditional_technique import *

if __name__ == '__main__':
	print('JACCARD ON K-MERS OR K-FINGERS')
	print('1: Calculate Jaccard similarity (exploring various k-mer size)')
	print('2: Calculate Jaccard similarity (fixed k-mer size)')
	print('3: Estimate Jaccard similarity (estimate on pairs of some files)')
	print('')
	print('MASH INTERFACING')
	print('4: Create new dataset from source')
	print('5: Preprocess all the files in a directory')
	print('6: Graph with MASH on k-mers and k-fingers (on preprocessed dataset) for short reads')
	print('7: Graph with MASH (and Jaccard) on k-mers and varying k-fingers (on preprocessed dataset) for long reads')
	print('')
	print('_: Quit')
	print('')
	action = input_int_not_neg('\nTell me what to do')

	if action == 1:
		seq1 = input_fasta_sequence('Tell me the first FASTA file')
		seq2 = input_fasta_sequence('Tell me the second FASTA file')

		print('Comparing:', seq1.get_name(), 'with', seq2.get_name())
		progressing_jaccard_on_kmers(seq1.get_data(), seq2.get_data(), 10, stdout_output)

	elif action == 2:
		files = input_files()
		kmer_size = input_kmer_size()

		sequences = [FastaSequence(file) for file in files]
		for seq1, seq2 in combinations(sequences, 2):
			monitor_jaccard_on_kmers_with_sequences(seq1, seq2, kmer_size, stdout_jaccard_on_kmers_output_function)

	elif action == 3:
		files = input_files()
		kmer_size = input_kmer_size()
		tolerance = input_float_from_zero_to_one('Tell me the tolerance', 0.15)

		great_estimation(files, kmer_size, tolerance)

	elif action == 4:
		source = input_fasta_sequence('Tell me the source for the creation of the new dataset (example: in/300.fasta)')
		dataset_count = input_int_not_neg('Tell me how many files do you want to create', 100)
		dataset_size_each_file = input_int_not_neg('Tell me the length of all the files', 600)
		kmer_size = input_kmer_size()
		factorization = input_factorization()
		kfinger_size = input_int_not_neg('Tell me the k-finger size', 3)
		use_super_fp = input_user_super_fp()
		min_overlap = input_int_not_neg('Tell me the minimum overlap size of two adjacent files', 250)
		upper_limit_for_fingers = input_int_not_neg('Tell me the max finger in the fingerprints of all files', 500)

		create_dataset_from_source(source, dataset_count, dataset_size_each_file, min_overlap, factorization,
		                           use_super_fp, upper_limit_for_fingers)

	elif action == 5:
		files = input_files_in_directory('Tell me the directory with the file to preprocess')
		factorization = input_factorization()
		use_super_fp = input_user_super_fp()
		split = input_int_not_neg('Tell me the split of the reads to apply (suggested 100 for long reads, 0 skips)', 0)

		preprocess_directory(files, factorization, use_super_fp, split)

	elif action == 6:
		files = input_files_in_directory('Tell me the dataset directory (that corresponds to the preprocessed dataset)')
		kmer_size = input_kmer_size()

		show_plot_mash_on_kmers_and_kfingers_based_on_preprocessed_dataset(files, kmer_size)

	elif action == 7:
		files = input_files_in_directory('Tell me the dataset directory (that corresponds to the preprocessed dataset)')
		kmer_size = input_kmer_size()
		factorization = input_factorization()
		kfinger_size = input_int_not_neg('Tell me the kfinger sizes to use on MASH', 4)
		use_super_fp = input_user_super_fp()
		split = input_int_not_neg('Tell me the split of the (long) reads to apply (0 skips the split)', 140)

		show_plot_mash_jaccard_varying_kfingers(files, kmer_size, factorization, kfinger_size, use_super_fp, split)

	else:
		print('Goodbye...️')
