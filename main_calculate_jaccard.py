from traditional_technique import jaccard_on_kmers
from Sequence import Sequence
from datetime import datetime


def progressing_jaccard_on_kmers(str1, str2, max_zero_similarities_before_stop, on_demand_output_function):
	kmer_size = 1
	p = jaccard_on_kmers(str1, str2, kmer_size)
	on_demand_output_function(kmer_size, p)
	kmer_size += 1
	last_zero_similarities_count = 0
	while last_zero_similarities_count < max_zero_similarities_before_stop:
		j = jaccard_on_kmers(str1, str2, kmer_size)
		if p == 0:
			last_zero_similarities_count += 1
		else:
			last_zero_similarities_count = 0
		p = j
		on_demand_output_function(kmer_size, j)
		kmer_size += 1


def calc_similarity_with_kmers_from_fasta_pair(seq1: Sequence, seq2: Sequence, k, *output_functions):
	before = datetime.now()
	similarity = jaccard_on_kmers(seq1.get_data(), seq2.get_data(), k)
	after = datetime.now()
	duration = after - before
	for output_function in output_functions:
		output_function(seq1.get_name(), seq2.get_name(), k, similarity, duration)


def stdout_output_function(name1, name2, k, similarity, duration):
	print("%2s <-> %2s  (k=%d) -> %6.2f" % (name1, name2, k, similarity * 100) + "%" + (
			"    (%d microseconds)" % duration.microseconds))
