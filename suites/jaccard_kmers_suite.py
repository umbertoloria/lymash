from files_manager import input_files
from traditional_technique import jaccard_on_kmers, monitor_jaccard_on_kmers_with_sequences, \
	stdout_monitor_jaccard_on_kmers_output_function


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


def stdout_output(kmer_size, similarity, prev_similarity=None):
	if prev_similarity is None:
		print("%3d => %10.9f" % (kmer_size, similarity))
	else:
		diff = (1 - similarity / prev_similarity) * 100
		print("%3d => %10.9f   (diminuisce del %5.3f" % (kmer_size, similarity, diff) + '%)')


def jaccard_kmers_suite_main(action):
	if action == 2:
		from sequences.Sequence import FastaSequence
		seq1 = FastaSequence(input("Tell me the first file> "))
		seq2 = FastaSequence(input("Tell me the second file> "))
		print('Comparing:', seq1.get_name(), 'with', seq2.get_name(), '\n')
		progressing_jaccard_on_kmers(seq1.get_data(), seq2.get_data(), 10, stdout_output)
	elif action == 3:
		from sequences.Sequence import FastaSequence
		from itertools import combinations
		files = input_files()
		kmer_size = input("Tell me the kmer size (usually 21)> ")
		kmer_size = 21 if kmer_size == "" else int(kmer_size)
		sequences = [FastaSequence(file) for file in files]
		for seq1, seq2 in combinations(sequences, 2):
			monitor_jaccard_on_kmers_with_sequences(seq1, seq2, kmer_size,
			                                        stdout_monitor_jaccard_on_kmers_output_function)
	else:
		return False
