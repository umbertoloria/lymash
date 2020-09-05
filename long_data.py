from Sequence import FastaSequence

import os

def create_long_data():
	dirname = "long_dataset"
	if not os.path.isdir(dirname):
		os.mkdir(dirname)

	for i in range(1,100):
		file = open(os.path.join(dirname + "/" + str(i) + ".fasta"), "w")
		file.write("@m #short read generated\n")
		for w in range(i*5,i*5+5):
			f = open("in/"+str(w)+".fasta", "r")
			lines = f.readlines()[1:]
			file.writelines(lines)


create_long_data()
