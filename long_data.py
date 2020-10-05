import os


def create_long_data():
	dirname = "long_dataset4"
	if not os.path.isdir(dirname):
		os.mkdir(dirname)

	for i in range(100,200):
		file = open(os.path.join(dirname + "/" + str(i) + ".fasta"), "w")
		file.write("@m54329U_190528_231241/17304174/ccs #short read generated\n")
		rand = 0
		for w in range(i*10, i*10+10):
			f = open("in/"+str(w)+".fasta", "r")
			lines = f.readlines()[1:]
			file.writelines(lines)


create_long_data()
