from factorization import subdivide


def get_reads_from_fq_file(filepath):
	f = open(filepath)
	i = 0
	line = True
	last_title = None
	while line:
		line = f.readline()
		if i % 4 == 0:
			last_title = line
			if last_title.endswith('\n'):
				last_title = last_title[:-1]
		if i % 4 == 1:
			yield last_title, line
		i += 1


def read_from_fasta(filepath):
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


def write_into_fasta(filepath, read, title=""):
	f = open(filepath, "w")
	f.write(">")
	f.write(title)
	f.write("\n")
	for line in subdivide(read, 70):
		f.write(line)
		f.write("\n")
	f.close()
