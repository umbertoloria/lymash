import os
from collections import defaultdict
from itertools import combinations

from new_combined_technique import get_csv_exporter, get_grafico_exporter, new_combined_technique_analyzer
from sequences import FastaSequence


def great_estimation(files: list, kmer_size: int, tolerance: float):
	"""Analyzes every pair of given files comparison calculation using Jaccard on k-fingers. Creates a .csv file for
	every pair of files, and a summary.txt file that summarizes the results of all comparisons. Also, for every pair of
	files, opens a comparison plot showing the comparisons survived to the 'tolerance' filter."""

	dirname = []
	for file in files:
		file = os.path.basename(file)
		if '.' in file:
			file = file[:file.rfind('.')]
		dirname.append(file)
	dirname = '-'.join(dirname) + ' kmer ' + str(kmer_size) + '-tolerance ' + str(int(tolerance * 100))
	attempt = 1
	if os.path.isdir(dirname):
		attempt = 2
		dirname = dirname + ' (2)'
	while os.path.isdir(dirname):
		attempt += 1
		dirname = dirname[:dirname.rfind('(')] + '(' + str(attempt) + ')'
	os.mkdir(dirname)

	sequences = [FastaSequence(file) for file in files]

	# Exporters
	csv_exporter = get_csv_exporter(dirname)
	grafico_exporter = get_grafico_exporter

	# Summary
	summary = open(dirname + '/summary.txt', 'w')

	# Sopravvissuti
	sopravvissuti = defaultdict(int)

	for seq1, seq2 in combinations(sequences, 2):
		print('Comparing', seq1.get_name(), '<->', seq2.get_name())

		result = new_combined_technique_analyzer(seq1, seq2, kmer_size, tolerance, csv_exporter, grafico_exporter)

		for factorization, configs in result['factorizations'].items():
			for split, window_size in configs:
				sopravvissuti[(split, window_size, factorization)] += 1

		summary.write(seq1.get_name() + ' - ' + seq2.get_name() + '\n')
		summary.write(seq1.get_title() + '\n')
		summary.write(seq2.get_title() + '\n\n')
		summary.write('k-mer size ->    ' + str(kmer_size) + '\n')
		summary.write('jaccard   ->    %5.2f' % (result['jaccard'] * 100) + '%\n')
		summary.write('tolerance ->    %5.2f' % (tolerance * 100) + '%\n\n')
		for factorization, configs in result['factorizations'].items():
			summary.write('  %13s  %3d\n' % (factorization, len(configs)))
		summary.write('\n----------------------------------------------------------------------\n\n')

	sopravvissuti = sorted(sopravvissuti.items(), key=lambda p: p[1], reverse=True)
	for (split, window_size, factorization), count in sopravvissuti:
		summary.write('%3d-%1d %15s   -> %3d\n' % (split, window_size, factorization, count))
	summary.close()
