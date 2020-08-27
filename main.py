#!/bin/env python3

import main_calculate_jaccard
import main_estimate_jaccard
import main_export_fq_in_fasta_files
import main_mash_interface
import main_process_files_in_directory


def main_bello():
	pass


if __name__ == "__main__":
	print("1: Export FASTA from FQ")
	print("2: Calculate Jaccard similarity (exploring various kmer size)")
	print("3: Calculate Jaccard similarity (fixed kmer size)")
	print("4: Estimate Jaccard similarity (estimate on pairs of some files)")
	print("5: Estimate Jaccard similarity (just two given files)")
	print("6: Preprocess all the files in a directory")
	print("7: Encode in ASCII the fingers of two FASTA files")
	print("_: Quit")
	action = int(input("\nTell me what to do > "))
	if action == 1:
		main_export_fq_in_fasta_files.main()
	elif action == 2:
		main_calculate_jaccard.explore_kmer_size()
	elif action == 3:
		main_calculate_jaccard.fixed_kmer_size()
	elif action == 4:

		files = ["in/1.fasta", "in/2.fasta", "in/3.fasta"]  # , "in/4.fasta", "in/5.fasta", "in/6.fasta"]
		configs = []
		for i in range(len(files) - 1):
			for j in range(i + 1, len(files)):
				configs.append((files[i], files[j]))

		from collections import defaultdict

		coppie_famose = defaultdict(int)
		for file1, file2 in configs:
			print("Comparing", file1, "-", file2)
			main_estimate_jaccard.estimate_jaccard_from_filepaths(file1, file2, 0.1, coppie_famose,
			                                                      main_estimate_jaccard.export_csv)
		coppie_famose_sorted = sorted(coppie_famose.items(), key=lambda p: p[1], reverse=True)
		for coppia_famosa in coppie_famose_sorted:
			split, window_size, alg = coppia_famosa[0]
			print("%3d-%1d %15s -> %3d" % (split, window_size, alg, coppia_famosa[1]))

	elif action == 5:
		from collections import defaultdict

		file1 = input("Tell me the first file> ")
		file2 = input("Tell me the second file> ")

		if file1 == "":
			file1 = "in/2.fasta"

		if file2 == "":
			file2 = "in/3.fasta"

		coppie_famose = defaultdict(int)
		main_estimate_jaccard.estimate_jaccard_from_filepaths(file1, file2, 0.1, coppie_famose,
		                                                      main_estimate_jaccard.export_csv)
		coppie_famose_sorted = sorted(coppie_famose.items(), key=lambda p: p[1], reverse=True)
		for coppia_famosa in coppie_famose_sorted:
			split, window_size, alg = coppia_famosa[0]
			print("%3d-%1d %15s -> %3d" % (split, window_size, alg, coppia_famosa[1]))

	elif action == 6:
		main_process_files_in_directory.main()
	elif action == 7:
		main_mash_interface.main()
	elif action == 8:
		main_bello()
	else:
		print("Goodbye...️")
