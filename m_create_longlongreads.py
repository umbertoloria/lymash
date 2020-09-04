import os

from format import read_from_fasta
from factorization import subdivide


def main(dir, maxlen):
	files = os.listdir(dir)
	files = [dir + "/" + file for file in files]

	cur_source_index = 0
	cur_source_data = None

	cur_dest_index = 0
	cur_dest_data = ""

	while cur_source_index < len(files):
		if cur_source_data is None:
			cur_source_data = read_from_fasta(files[cur_source_index])[1]
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
