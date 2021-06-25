#!/usr/bin/env python

import sys
import gzip

from pathlib import Path

import pybedtools 


def main():
	if len(sys.argv) < 7:
		print ("Usage : " + sys.argv[0] + " [chrominfo_filename] [bin_size] [slop_size] [gap_bedgz_filename] [blacklist_bedgz_filename] [filtered_genomic_bin_bedgz_filename]")
		sys.exit()
	else:
		chrominfo_filename = sys.argv[1]
		bin_size = int(sys.argv[2])
		slop_size = int(sys.argv[3])
		gap_bedgz_filename = sys.argv[4]
		blacklist_bedgz_filename = sys.argv[5]
		filtered_genomic_bin_bedgz_filename = sys.argv[6]
		make_filtered_genomic_bin(chrominfo_filename,bin_size,slop_size,gap_bedgz_filename , blacklist_bedgz_filename , filtered_genomic_bin_bedgz_filename)


def make_filtered_genomic_bin(chrominfo_filename,bin_size,slop_size,gap_bedgz_filename , blacklist_bedgz_filename , filtered_genomic_bin_bedgz_filename):
	genomic_bin_bed = pybedtools.BedTool().window_maker(g=chrominfo_filename,w=bin_size+slop_size*2,s=bin_size,stream=True)

	gap_bed = pybedtools.BedTool(gap_bedgz_filename)
	blacklist_bed = pybedtools.BedTool(blacklist_bedgz_filename)
	filter_bed = gap_bed.cat(blacklist_bed).filter(lambda interval : not( interval.chrom in {"chrY","chrM"}) )
	filtered_genomic_bin_bed = genomic_bin_bed.intersect(filter_bed,v=True)
	filtered_genomic_bin_bed = filtered_genomic_bin_bed.filter(lambda interval : len(interval) == bin_size+slop_size*2)
	Path.mkdir(Path(filtered_genomic_bin_bedgz_filename).parent,parents=True,exist_ok=True) 
	write_genomic_bin(filtered_genomic_bin_bed,filtered_genomic_bin_bedgz_filename)


def write_genomic_bin(filtered_genomic_bin_bed,filtered_genomic_bin_bedgz_filename):
	with gzip.open(filtered_genomic_bin_bedgz_filename,"wt") as filtered_genomic_bin_bedgz_file:
		for index,interval in enumerate(filtered_genomic_bin_bed):
			interval_name = "bin_" + str(index+1).zfill(6) 
			filtered_genomic_bin_bedgz_file.write("\t".join(map(str,list(interval))) + "\t" + interval_name + "\n" )


if __name__ == "__main__":
	main()
