from traditional_technique import jaccard_on_kmers


def progressing_jaccard_on_kmers(str1, str2, max_zero_similarities_before_stop, on_demand_output_function):
	kmer_size = 1
	p = jaccard_on_kmers(str1, str2, kmer_size)
	on_demand_output_function(kmer_size, p)
	kmer_size += 1
	last_zero_similarities_count = 0
	while last_zero_similarities_count < max_zero_similarities_before_stop:
		j = jaccard_on_kmers(str1, str2, kmer_size)
		if p == 0:
			last_zero_similarities_count += 1
		else:
			last_zero_similarities_count = 0
		on_demand_output_function(kmer_size, j, p)
		p = j
		kmer_size += 1


def stdout_output(kmer_size, similarity, prev_similarity=None):
	if prev_similarity is None:
		print('%3d => %10.9f' % (kmer_size, similarity))
	else:
		if prev_similarity > 0:
			diff = (1 - similarity / prev_similarity) * 100
		else:
			diff = similarity * 100
		print('%3d => %10.9f   (diminuisce del %5.3f' % (kmer_size, similarity, diff) + '%)')
