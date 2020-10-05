from sliding_windows import ksliding
from jaccard import jaccard
from sequences.Sequence import Sequence
from datetime import datetime


def jaccard_on_kmers(str1, str2, k):
	return jaccard(set(ksliding(str1, k)), set(ksliding(str2, k)))


def monitor_jaccard_on_kmers_with_sequences(seq1: Sequence, seq2: Sequence, k, *output_functions):
	before = datetime.now()
	similarity = jaccard_on_kmers(seq1.get_data(), seq2.get_data(), k)
	after = datetime.now()
	duration = after - before
	for output_function in output_functions:
		output_function(seq1.get_name(), seq2.get_name(), k, similarity, duration)


def stdout_monitor_jaccard_on_kmers_output_function(name1, name2, k, similarity, duration):
	print("%2s <-> %2s  (k=%d) -> %6.2f" % (name1, name2, k, similarity * 100) + "%" + (
			"    (%d microseconds)" % duration.microseconds))
