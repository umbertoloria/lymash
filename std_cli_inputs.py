from cli_input import input_int_not_neg, input_file, input_bool
from sequences import FastaSequence


def input_kmer_size(message: str = 'Tell me the k-mer size'):
	return input_int_not_neg(message, 21)


def input_fasta_sequence(message: str):
	return FastaSequence(input_file(message, '.fasta'))


def input_user_super_fp():
	return input_bool('Tell me if you want to use the super-fingerprint')
