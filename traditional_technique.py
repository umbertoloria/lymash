from sliding_windows import ksliding
from jaccard import jaccard


def jaccard_on_kmers(str1, str2, k):
	return jaccard(set(ksliding(str1, k)), set(ksliding(str2, k)))
