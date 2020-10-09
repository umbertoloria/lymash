#!/bin/env python3

import os
import main_export_fq_in_fasta_files
from sequences.Sequence import FastaSequence
from suites import jaccard_kmers_suite, mash_interfacing_suite, jaccard_kfingers_suite

if __name__ == "__main__":
	print("FORMATS")
	print("1: Export FASTA from FQ")
	print("")
	print("CALCULATE SIMILARITY: JACCARD ON K-MERS")
	print("2: Calculate Jaccard similarity (exploring various kmer size)")
	print("3: Calculate Jaccard similarity (fixed kmer size)")
	print("")
	print("ESTIMATE SIMILARITY: JACCARD ON K-FINGERS")
	print("4: Estimate Jaccard similarity (estimate on pairs of some files)")
	print("")
	print("MASH INTERFACING")
	print("5: Preprocess all the files in a directory (only short reads)")  # TODO:only short_reads???
	print("7: Graph with MASH on k-mers and k-fingers (on preprocessed dataset)")
	print("8: Graph with MASH and only-Jaccard on k-mers and varying k-fingers")
	print("9: Crea dataset di long long reads")
	print("")
	print("BOH...")
	print("_: Quit")
	print("")
	action = int(input("\nTell me what to do > "))

	if jaccard_kmers_suite.jaccard_kmers_suite_main(action) != False:
		pass
	elif jaccard_kfingers_suite.jaccard_kfingers_suite_main(action) != False:
		pass
	elif mash_interfacing_suite.mash_interfacing_suite_main(action) != False:
		pass

	elif action == 1:
		main_export_fq_in_fasta_files.main()

	elif action == 9:

		# TODO: testa bene questa cosa!
		def generate_parts(dir, maxlen):
			files = os.listdir(dir)
			files = [dir + "/" + file for file in files]

			cur_source_index = 0
			cur_source_data = None

			cur_dest_index = 0
			cur_dest_data = ""

			while cur_source_index < len(files):
				if cur_source_data is None:
					seq = FastaSequence(files[cur_source_index])
					cur_source_data = seq.get_data()
				data_to_consume_len = len(cur_source_data)
				if len(cur_dest_data) + data_to_consume_len <= maxlen:
					cur_dest_data += cur_source_data
					cur_source_data = None
					cur_source_index += 1
				else:
					space_to_fill = maxlen - len(cur_dest_data)
					cur_dest_data += cur_source_data[:space_to_fill]
					cur_source_data = cur_source_data[space_to_fill:]
					cur_dest_index += 1
					yield cur_dest_data
					cur_dest_data = ""
			if len(cur_dest_data) > 0:
				yield cur_dest_data


		dir = "prova"
		maxlen = 250
		for part in generate_parts(dir, maxlen):
			print(part)

	else:
		print("Goodbye...️")
# TODO: cancellare cartella DSA? E quali altre?
