from format import read_from_fasta


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
		self.__title, self.__data = read_from_fasta(fasta_path)

	def get_name(self) -> str:
		return self.__name

	def get_title(self) -> str:
		return self.__title

	def get_data(self) -> str:
		return self.__data
