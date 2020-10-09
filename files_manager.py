import os


# TODO: creare una bella libreria per gestire tutti gli input (anche l'input di UNA SOLA directory)
def input_files():
	print("Tell me files to compare ")
	files = []
	path = input("> ")
	while path != "":
		if not os.path.exists(path):
			print("File o directory non esistente")
		elif os.path.isdir(path):
			new_files = os.listdir(path)
			new_files.sort()
			for new_file in new_files:
				new_path = path + "/" + new_file
				if new_path not in files:
					files.append(new_path)
		elif path not in files:
			files.append(path)
		path = input("> ")
	return files
