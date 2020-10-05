def _read_from_fasta(filepath):
	f = open(filepath, "r")
	line = f.readline()
	if line.endswith('\n'):
		line = line[:-1]
	title = line[1:]
	read = ""
	while line:
		line = f.readline()
		if line.endswith('\n'):
			line = line[:-1]
		read += line
	return title, read


class Sequence:

	def get_name(self):
		raise NotImplemented()

	def get_title(self):
		raise NotImplemented()

	def get_data(self):
		raise NotImplemented()


class FastaSequence(Sequence):

	def __init__(self, fasta_path):
		self.__name = fasta_path.split("/")[-1]
		self.__title, self.__data = _read_from_fasta(fasta_path)

	def get_name(self) -> str:
		return self.__name

	def get_title(self) -> str:
		return self.__title

	def get_data(self) -> str:
		return self.__data
