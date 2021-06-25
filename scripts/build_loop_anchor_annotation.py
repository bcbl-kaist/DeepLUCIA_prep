#!/usr/bin/env python

import sys
import gzip
import itertools 
import collections

from pathlib import Path

def main():
	if len(sys.argv) < 5:
		print ("Usage : " + sys.argv[0] + " [filtered_genomic_bin_bedgz_filename] [loop_bedpegz_dirname] [bin_size] [anchor_annotation_txtgz_filename] [loop_annotation_txtgz_filename]")
		sys.exit()
	else:
		filtered_genomic_bin_bedgz_filename = sys.argv[1]
		loop_bedpegz_dirname = sys.argv[2]
		bin_size = int(sys.argv[3])
		anchor_annotation_txtgz_filename = sys.argv[4]
		loop_annotation_txtgz_filename = sys.argv[5]
		build_loop_anchor_annotation(filtered_genomic_bin_bedgz_filename, loop_bedpegz_dirname , bin_size  , anchor_annotation_txtgz_filename , loop_annotation_txtgz_filename )

def build_loop_anchor_annotation(filtered_genomic_bin_bedgz_filename, loop_bedpegz_dirname , bin_size  , anchor_annotation_txtgz_filename , loop_annotation_txtgz_filename ):
	bin_id_to_index,index_to_bin_id = load_genomic_bin(filtered_genomic_bin_bedgz_filename)
	sample_to_index_pair_set = load_loop_bedpe(bin_id_to_index,loop_bedpegz_dirname , bin_size)
	write_anchor_annotation(index_to_bin_id,sample_to_index_pair_set,anchor_annotation_txtgz_filename)
	write_loop_annotation(index_to_bin_id,sample_to_index_pair_set,loop_annotation_txtgz_filename)


def write_loop_annotation(index_to_bin_id,sample_to_index_pair_set,loop_annotation_txtgz_filename):
	Path.mkdir(Path(loop_annotation_txtgz_filename).parent, parents=True, exist_ok=True)
	index_pair_to_sample_set = collections.defaultdict(set)
	for sample,index_pair_set in sample_to_index_pair_set.items():
		for index_pair in index_pair_set:
			index_pair_to_sample_set[index_pair].add(sample)

	with gzip.open(loop_annotation_txtgz_filename,"wt") as loop_annotation_txtgz_file:
		loop_annotation_txtgz_file.write("anchor_id_one\tanchor_id_two\tanchor_index_one\tanchor_index_two\tcorresponding_sample_id_label\n")
		for index_pair in sorted(index_pair_to_sample_set.keys()):
			index_one,index_two = index_pair
			bin_id_one = index_to_bin_id[index_one] 
			bin_id_two = index_to_bin_id[index_two]
			sample_label = ",".join(sorted(index_pair_to_sample_set[index_pair]))
			loop_annotation_txtgz_file.write("\t".join(map(str,[bin_id_one,bin_id_two,index_one,index_two,sample_label])) + "\n")


def write_anchor_annotation(index_to_bin_id,sample_to_index_pair_set,anchor_annotation_txtgz_filename):
	Path.mkdir(Path(anchor_annotation_txtgz_filename).parent, parents=True, exist_ok=True)
	index_to_sample_set = collections.defaultdict(set)
	for sample,index_pair_set in sample_to_index_pair_set.items():
		for index_pair in index_pair_set:
			index_one,index_two = index_pair
			index_to_sample_set[index_one].add(sample)
			index_to_sample_set[index_two].add(sample)

	with gzip.open(anchor_annotation_txtgz_filename,"wt") as anchor_annotation_txtgz_file:
		anchor_annotation_txtgz_file.write("anchor_index\tanchor_id\tcorresponding_sample_id_label\n")
		for index in sorted(index_to_bin_id.keys()):
			bin_id = index_to_bin_id[index]
			sample_label = ",".join(sorted(index_to_sample_set[index])) if index in index_to_sample_set else "None"
			anchor_annotation_txtgz_file.write("\t".join(map(str,[index,bin_id,sample_label])) + "\n")


def load_loop_bedpegz_dir(loop_bedpegz_dirname):
	sample_to_loop_bedpegz_filename = {}
	for loop_bedpegz_filename in Path(loop_bedpegz_dirname).glob("*/*.bedpe.gz"):
		sample = loop_bedpegz_filename.name.split("_")[0]
		sample_to_loop_bedpegz_filename[sample] = loop_bedpegz_filename

	return sample_to_loop_bedpegz_filename



def load_loop_bedpe(bin_id_to_index,loop_bedpegz_dirname , bin_size):
	sample_to_index_pair_set = {}
	sample_to_loop_bedpegz_filename = load_loop_bedpegz_dir(loop_bedpegz_dirname)
	for sample,loop_bedpegz_filename in sample_to_loop_bedpegz_filename.items():
		sample_to_index_pair_set[sample] = set([])
		with gzip.open(loop_bedpegz_filename,"rt") as loop_bedpegz_file:
			for rawline in itertools.islice(loop_bedpegz_file,1,None):
				fields = rawline.strip().split()
				if fields[0] == fields[3]:
					chrom = "chr" + fields[0]
					bin_one_start_gen = range(int(fields[1]),int(fields[2]),bin_size)
					bin_two_start_gen = range(int(fields[4]),int(fields[5]),bin_size)
					for bin_one_start,bin_two_start in itertools.product(bin_one_start_gen,bin_two_start_gen):
						bin_id_one = chrom + ":" + str(bin_one_start) + "-" + str(bin_one_start+bin_size)
						bin_id_two = chrom + ":" + str(bin_two_start) + "-" + str(bin_two_start+bin_size)
						#print(bin_id_one , bin_id_one in bin_id_to_index , bin_id_two , bin_id_two in bin_id_to_index)

						if bin_id_one in bin_id_to_index and bin_id_two in bin_id_to_index:
							index_one = bin_id_to_index[bin_id_one]
							index_two = bin_id_to_index[bin_id_two]
							index_pair = (index_one,index_two)
							sample_to_index_pair_set[sample].add(index_pair)

	return sample_to_index_pair_set


def load_genomic_bin(filtered_genomic_bin_bedgz_filename):
	bin_id_to_index = {}
	index_to_bin_id = {}

	with gzip.open(filtered_genomic_bin_bedgz_filename,"rt") as filtered_genomic_bin_bedgz_file:
		for rawline in filtered_genomic_bin_bedgz_file:
			fields = rawline.strip().split()
			bin_id = fields[0] + ":" + fields[1] + "-" + fields[2]
			index = int(fields[3].split("_")[1])
			bin_id_to_index[bin_id] = index
			index_to_bin_id[index] = bin_id

	return bin_id_to_index,index_to_bin_id


if __name__ == "__main__":
	main()



