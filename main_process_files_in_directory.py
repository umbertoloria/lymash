import os

from factorization import get_fingers_after_text_subdividing, subdivide
from format import read_from_fasta


def main():
	dir = input("Tell me the directory with the file to preprocess> ")
	factorization_method = input(
		"Tell me the algorithm for factoring\n   (cfl,icfl,cfl_icfl,cfl_comb,icfl_comb,cfl_icfl_comb)> ")
	subdivision = int(input("Tell me the subdivision lenght (max 200 for ASCII sake)> "))

	if not os.path.isdir(dir + "_preprocessed"):
		os.mkdir(dir + "_preprocessed")
	for root, dirs, filepaths in os.walk(dir):
		for filepath in filepaths:
			if filepath.endswith(".fasta"):
				title, read = read_from_fasta(dir + "/" + filepath)
				fingers = get_fingers_after_text_subdividing(factorization_method, read, subdivision)
				file = open(os.path.join(dir + "_preprocessed", "preprocessed_" + filepath), "w")
				file.write(title + " preprocessed encoding lyndon fact lengths with ascii code\n")
				for fingers_row in subdivide(fingers, 70):
					for finger in fingers_row:
						num = ord('0') - 1 + finger
						if num > ord('Z'):
							print(" *** using not printable character (" + str(num) + "):", chr(num))
						file.write(chr(num))
					file.write("\n")
				file.close()
