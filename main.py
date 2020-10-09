#!/bin/env python3

from suites import jaccard_kmers_suite, mash_interfacing_suite, jaccard_kfingers_suite

if __name__ == "__main__":
	print("CALCULATE SIMILARITY: JACCARD ON K-MERS")
	print("1: Calculate Jaccard similarity (exploring various kmer size)")
	print("2: Calculate Jaccard similarity (fixed kmer size)")
	print("")
	print("ESTIMATE SIMILARITY: JACCARD ON K-FINGERS")
	print("3: Estimate Jaccard similarity (estimate on pairs of some files)")
	print("")
	print("MASH INTERFACING")
	print("4: Preprocess all the files in a directory (only short reads)")
	print("5: Graph with MASH on k-mers and k-fingers (on preprocessed dataset) for short reads")
	print("6: Graph with MASH (and Jaccard) on k-mers and varying k-fingers (on preprocessed dataset) for long reads")
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
	else:
		print("Goodbye...️")
