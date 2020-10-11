class Sequence:

	def get_name(self):
		raise NotImplemented()

	def get_title(self):
		raise NotImplemented()

	def get_data(self):
		raise NotImplemented()


class FastaSequence(Sequence):

	def __init__(self, fasta_path):
		self.__name = fasta_path.split('/')[-1]
		f = open(fasta_path, 'r')
		line = f.readline()
		if line.startswith('>'):
			line = line[1:]
		if line.endswith('\n'):
			line = line[:-1]
		self.__title = line
		self.__data = ''
		line = True
		while line:
			line = f.readline()
			if line.endswith('\n'):
				line = line[:-1]
			self.__data += line

	def get_name(self) -> str:
		return self.__name

	def get_title(self) -> str:
		return self.__title

	def get_data(self) -> str:
		return self.__data
