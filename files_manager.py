import os


def input_files_in_directory(message: str):
	directory = input(message + '> ')
	while not os.path.isdir(directory):
		print('   Directory not found')
		directory = input('> ')
	files = [directory + '/' + file for file in os.listdir(directory)]
	files.sort()
	return files


def input_factorization(facts=None):
	if facts is None:
		facts = ['cfl', 'icfl', 'cfl_icfl', 'cfl_comb', 'icfl_comb', 'cfl_icfl_comb']
	print('Tell me the factorization method')
	factorization = input('   (' + ','.join(facts) + ')> ')
	factorization = factorization.strip().lower()
	while factorization not in facts:
		print('   Factorization method not suitable')
		factorization = input('> ')
		factorization = factorization.strip().lower()
	return factorization


def input_files():
	print('Tell me files to compare ')
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
