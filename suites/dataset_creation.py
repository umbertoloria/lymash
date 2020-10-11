import os
import random

from new_combined_technique import create_fingerprint, subdivide
from sequences import Sequence


def create_dataset_from_source(source: Sequence, count: int, size: int, min_overlap: int, factorization: str,
                               use_super_fp: bool, upper_limit_for_fingers: int):
	""" This function creates a dataset with specific properties. The number of files in it is defined by the 'count'
	parameter, and the character count by the 'size' parameter. Two adjacents files (for example 1.fasta and 2.fasta)
	shares an overlap, and all the overlaps are randomly defined from a fixed number of characters up to all. This
	minimum fixed number of characters for every overlap is defined by the 'min_overlap' parameter. The most important
	property of the generated dataset is that the fingerprints (whose model depends on 'factorization' and
	'use_super_fp' parameters) contains fingers less or equal then the 'upper_limit_for_fingers' parameter."""

	dirname = 'dataset_' + source.get_name().split('.')[0]
	if not os.path.isdir(dirname):
		os.mkdir(dirname)

	dataset = ['' for _ in range(count)]

	total_max = upper_limit_for_fingers
	while total_max >= upper_limit_for_fingers:
		start = 0
		for i in range(count):
			dataset[i] = source.get_data()[start:start + size]
			rand = random.randrange(min_overlap, size)
			start += size - rand
		total_max = 0
		for text in dataset:
			fingerprint = create_fingerprint(text, factorization, use_super_fp)
			cur_max = max(fingerprint)
			if cur_max > total_max:
				total_max = cur_max
		print('Max finger: {} '.format(total_max))

	for i in range(count):
		file = open(dirname + '/' + str(i + 1) + '.fasta', 'w')
		file.write('>\'' + source.get_title() + '\' short read randomized part\n')
		text = dataset[i]
		for part in subdivide(text, 70):
			file.write(part)
			file.write('\n')
		file.close()
