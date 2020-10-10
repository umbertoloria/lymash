import os
import random

from files_manager import input_factorization
from new_combined_technique import create_fingerprint, subdivide
from sequences.Sequence import FastaSequence

# TODO ctrl-o non solo ctrl-l!
# TODO: another single main


source_filepath = input('Tell me the source from which extrapolate parts to form the new dataset (in/300.fasta)> ')
if source_filepath == '':
	source_filepath = 'in/300.fasta'
source = FastaSequence(source_filepath)

dataset_count = input('Tell me how many files do you want the new dataset to contain (100)> ')
dataset_count = 100 if dataset_count == '' else int(dataset_count)

dataset_size_each_file = input('Tell me the length of all files in the new dataset (600)> ')
dataset_size_each_file = 600 if dataset_size_each_file == '' else int(dataset_size_each_file)

# TODO: automate input kmer/kfingers
kmer_size = input('Tell me the k-mer size (usually 21)> ')
kmer_size = 21 if kmer_size == '' else int(kmer_size)

factorization = input_factorization()

kfinger_size = input('Tell me the k-finger size (3)> ')
kfinger_size = 3 if kfinger_size == '' else int(kfinger_size)

use_super_fp = input('Tell me if you want to use the super-fingerprint [Y/n]> ')
use_super_fp = use_super_fp == 'Y' or use_super_fp == 'y' or use_super_fp == ''

min_overlap = input('Tell me the minimum overlap size two adjacent files will have in this new dataset (250)> ')
min_overlap = 250 if min_overlap == '' else int(min_overlap)

print()

dirname = 'dataset_' + source.get_name().split('.')[0]
if not os.path.isdir(dirname):
	os.mkdir(dirname)

dataset = ['' for i in range(dataset_count)]

total_max = 500
while total_max >= 500:
	start = 0
	for i in range(dataset_count):
		dataset[i] = source.get_data()[start:start + dataset_size_each_file]
		rand = random.randrange(min_overlap, dataset_size_each_file)
		start += dataset_size_each_file - rand

	total_max = 0
	for text in dataset:
		fingerprint = create_fingerprint(text, factorization, use_super_fp)
		cur_max = max(fingerprint)
		if cur_max > total_max:
			total_max = cur_max

	print('Max finger: {} '.format(total_max))

for i in range(dataset_count):
	file = open(dirname + '/' + str(i + 1) + '.fasta', 'w')
	file.write('>\'' + source.get_title() + '\' short read randomized part\n')
	text = dataset[i]
	for part in subdivide(text, 70):
		file.write(part)
		file.write('\n')
	file.close()
