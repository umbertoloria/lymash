from sequences.Sequence import Sequence
from sliding_windows import ksliding
from jaccard import jaccard
from traditional_technique import jaccard_on_kmers
from factorization import get_fingers_after_split


def estimate_jaccard(seq1: Sequence, seq2: Sequence, factorization: str, kfinger_size: int, use_super_fp: bool = False):
	from factorization import get_factors
	from lyndon.utils import _complement

	# Create (super-)fingerprints
	afp = [len(f) for f in get_factors(factorization, seq1.get_data())]
	if use_super_fp:
		afp += [0]
		afp += [len(f) for f in get_factors(factorization, _complement(seq1.get_data()))]
	afp = tuple(afp)

	bfp = [len(f) for f in get_factors(factorization, seq2.get_data())]
	if use_super_fp:
		bfp += [0]
		bfp += [len(f) for f in get_factors(factorization, _complement(seq2.get_data()))]
	bfp = tuple(bfp)

	# Getting k-fingers
	a_kfingers = set(ksliding(afp, kfinger_size, False))
	b_kfingers = set(ksliding(bfp, kfinger_size, False))

	jaccard_estimated = jaccard(a_kfingers, b_kfingers)
	return jaccard_estimated


def estimate_jaccard_difference_split(seq1: Sequence, seq2: Sequence, kmer_size: int, factorization, kfinger_size: int):
	from lyndon.utils import _complement

	read1 = seq1.get_data()
	read2 = seq2.get_data()

	jaccard_calculated = jaccard_on_kmers(read1, read2, kmer_size)

	split = 140

	a_factors_lengths = tuple(get_fingers_after_split(factorization, read1, split)) + (0,) + tuple(
		get_fingers_after_split(factorization, _complement(read1), split))
	b_factors_lengths = tuple(get_fingers_after_split(factorization, read2, split)) + (0,) + tuple(
		get_fingers_after_split(factorization, _complement(read2), split))

	a_fingerprint = set(ksliding(a_factors_lengths, kfinger_size, False))
	b_fingerprint = set(ksliding(b_factors_lengths, kfinger_size, False))

	jaccard_estimated = jaccard(a_fingerprint, b_fingerprint, True)

	return jaccard_calculated, jaccard_estimated
