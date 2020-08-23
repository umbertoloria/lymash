from factorization import subdivide
from format import get_reads_from_fq_file


def main():
	filepath = input("Tell me the filepath of the FQ file to \"explode\" > ")
	reads = get_reads_from_fq_file(filepath)
	i = 1
	for title, read in reads:
		f = open("in/" + str(i) + ".fasta", "w")
		i += 1
		parts = subdivide(read, 70)
		f.write(">")
		f.write(title)
		for part in parts:
			f.write("\n")
			f.write(part)
		f.close()
