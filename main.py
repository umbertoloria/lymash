#!/bin/env python3

import main_calculate_jaccard
import main_estimate_jaccard
import main_export_fq_in_fasta_files
import main_mash_interface
import main_process_files_in_directory

if __name__ == "__main__":
	print("1: Export FASTA from FQ")
	print("2: Calculate Jaccard similarity")
	print("3: Estimate Jaccard similarity")
	print("4: Preprocess all the files in a directory")
	print("5: Encode in ASCII the fingers of two FASTA files")
	print("_: Quit")
	action = int(input("\nTell me what to do > "))
	if action == 1:
		main_export_fq_in_fasta_files.main()
	elif action == 2:
		main_calculate_jaccard.main()
	elif action == 3:
		main_estimate_jaccard.main()
	elif action == 4:
		main_process_files_in_directory.main()
	elif action == 5:
		main_mash_interface.main()
	else:
		print("Goodbye...️")
