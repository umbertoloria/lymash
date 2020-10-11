import os

from sequences import FastaSequence


def input_kmer_size(message: str = 'Tell me the k-mer size'):
	return input_int_not_neg(message, 21)


def input_fasta_sequence(message: str):
	val = input(message + '> ')
	while not os.path.isfile(val) or not val.endswith('.fasta'):
		val = input('   invalid file> ')
	return FastaSequence(val)


def input_user_super_fp():
	val = input('Tell me if you want to use the super-fingerprint [Y/n]> ').lower()
	while val != 'y' and val != '' and val != 'n':
		val = input('   enter \'y\' (or nothing) for yes, \'n\' for no> ')
	return val != 'n'


def input_float_from_zero_to_one(message: str, default: float = None):
	num = input(message + (' (' + str(default) + ')' if default is not None else '') + '> ')
	while True:
		if num == '' and default is not None:
			return default
		elif num.isdecimal():
			return float(num)
		num = input('   enter value in [0, 1]> ')


def input_int_not_neg(message: str, default: int = None):
	num = input(message + (' (' + str(default) + ')' if default is not None else '') + '> ')
	while True:
		if num == '' and default is not None:
			return default
		elif num.isdigit():
			return int(num)
		num = input('   enter some digits> ')


def input_files_in_directory(message: str):
	directory = input(message + '> ')
	while not os.path.isdir(directory):
		print('   Directory not found')
		directory = input('> ')
	files = [directory + '/' + file for file in os.listdir(directory)]
	files.sort()
	return files


def input_factorization(message: str = 'Tell me the factorization method', facts=None, default: str = 'cfl_icfl_comb'):
	if facts is None:
		facts = ['cfl', 'icfl', 'cfl_icfl', 'cfl_comb', 'icfl_comb', 'cfl_icfl_comb']
	print(message + (' (' + default + ')' if default in facts else ''))
	factorization = input('   (' + ','.join(facts) + ')> ')
	factorization = factorization.strip().lower()
	while factorization not in facts:
		if default in facts and factorization == '':
			return default
		print('   Factorization method not suitable')
		factorization = input('> ')
		factorization = factorization.strip().lower()
	return factorization


def input_files(message: str = 'Tell me files to compare '):
	print(message, end='')
	files = []
	path = input('> ')
	while path != '':
		if not os.path.exists(path):
			print('   File o directory non esistente')
		elif os.path.isdir(path):
			new_files = os.listdir(path)
			new_files.sort()
			for new_file in new_files:
				new_path = path + '/' + new_file
				if new_path not in files:
					files.append(new_path)
		elif path not in files:
			files.append(path)
		path = input('> ')
	return files
