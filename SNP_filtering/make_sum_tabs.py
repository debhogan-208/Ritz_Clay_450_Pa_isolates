#!/usr/bin/env python3

# AUTHOR		:	ALAN COLLINS
# VERSION		:	v1
# DATE			:	2021-1-5
# DESCRIPTION	:	Make summary tables describing coverage and consensus score for each SNP of each isolate.

import sys
import os
import argparse


parser = argparse.ArgumentParser(
	description="Make summary tables describing coverage and consensus score for each SNP of each isolate in the input directory.")
parser.add_argument(
	"-i", dest="indir", required = True,
	help="Path to input directory containing snp_support files output by get_read_info_at_pos.py."
	)
parser.add_argument(
	"-o", dest="outfile", required = True,
	help="Path to output file prefix. _coverage.txt and _consensus_score.txt will be added to whatever you provide."
	)

args = parser.parse_args(sys.argv[1:])

indir = args.indir + '/' if args.indir[-1] != '/' else args.indir

snp_list = ['']
cov_to_write = []
con_score_to_write = []

for i, file in enumerate(os.listdir(indir)):
	infile = indir + file
	cov_list = [file[:-16]] # file path up to _snp_support.txt
	cons_list = [file[:-16]] # file path up to _snp_support.txt
	with open(infile, 'r') as fin:
		for j, line in enumerate(fin.readlines()[1:]):
			bits = line.split()
			if i == 0:
				snp_list.append(bits[2])
			if bits[2] == snp_list[j+1]:
				cov_list.append(bits[8])
				cons_list.append(bits[9])
			else:
				print("different SNP order in {} for SNP {}, {}".format(infile, snp_list[j], bits[2]))
				break
	if i == 0:
		cov_to_write.append(snp_list)
		con_score_to_write.append(snp_list)
	cov_to_write.append(cov_list)
	con_score_to_write.append(cons_list)

with open(args.outfile + "_coverage.txt", 'w+') as fout:
	fout.write("\n".join(["\t".join(line) for line in cov_to_write]) + '\n')

with open(args.outfile + "_consensus_score.txt", 'w+') as fout:
	fout.write("\n".join(["\t".join(line) for line in con_score_to_write]) + '\n')

