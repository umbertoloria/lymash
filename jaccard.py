def jaccard(a, b, verboose=False):
	intersection = len(a.intersection(b))
	union = len(a.union(b))
	if verboose:
		print("Jaccard: ", intersection, '/', union)
	return intersection / union
