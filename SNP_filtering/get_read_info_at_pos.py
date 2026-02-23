#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v1
# DATE        :  2021-1-8
# DESCRIPTION :  Look up SNP locations in reads in .sam files and report information about the reads at those locations in a file with the following columns:
# 1.Contig	2.Position	3.Contig_Position	4.Reference_Reads_agree	5.Assembly_Reads_agree	6.Reference_base	7.Isolate_base	8.Assembly_base	9.Coverage	10.Consensus_score	11.Read_base_calls	12.Read_qualities	13.Read_identifiers
# 1. The contig in which the SNP is found.
# 2. The position within that contig where the SNP is found (base 1 position - i.e. starts counting at position 1).
# 3. The contig number and position in format Contig_Position.
# 4. Binary indicator of whether the reference genome and the majority base call in the reads is the same. 0 = different, 1 = same.
# 5. Binary indicator of whether the base called by breseq and the majority base call in the reads is the same. 0 = different, 1 = same.
# 6. The base in the reference genome at the position.
# 7. The majority base call in the reads at the position.
# 8. The base called by breseq at the position (if known).
# 9. The number of reads that align at that position (removes reads that align nearby but are clipped according to CIGAR string in .sam file.
# 10. Consensus score = (sum of quality scores for base calls agreeing with majority call) / (sum of all quality scores at this position among all reads).
# 11. Ordered list of the bases found in the reads aligned at SNP position. Order of this list corresponds to the order of the qualities and identifiers lists in columns 12 and 13.
# 12. Ordered list of the quality scores of base calls in reads aligned at SNP position. Quality scores are converted from phred-33 to integer form for human-readability.
# 13. Ordered list of read names/identifiers for reads that align to SNP position


import argparse
import sys
import re
import time
import multiprocessing
import math
import os
import gzip
import pickle
import snpclasses

start_time = time.time()

def fasta_to_dict(FASTA_file):
	""" Convert fasta sequence into python dict object of form {fasta_header: sequence}
	
		Args:
			FASTA_file (str): fasta format DNA sequence
	
		Returns:
			dict: dict with fasta headers as keys and associated seuqences as values
	"""
	fasta_dict = {}
	fastas = FASTA_file.split(">")
	trimmed_fastas = []
	for i in fastas:
		if len(i) != 0:
			trimmed_fastas.append(i)

	fastas = trimmed_fastas

	for i in fastas:
		header = i.split("\n")[0]
		seq = "".join(i.split("\n")[1:])
		fasta_dict[header] = seq

	return fasta_dict


def cigar_mod_read(seq1, qual1, cigar):
	""" According to the M, I, and D components of a cigar string, modifies a seq and quality score strings so that they are in register with reference sequence 
	
		Args:
			seq1 (str): nucleotide sequence of read.
			qual1 (str): quality score string of read.
			cigar (str): cigar string of sam alignment of read against reference genome.
	
		Returns:
			tuple: seq2 (str) cigar-modified sequence of read,
				   qual2 (str) cigar-modified quality scores
	"""
	cig_sects = re.findall('([0-9]+[A-Z]+)', cigar)

	seq2 =''
	qual2 = ''

	count = 0

	for sect in cig_sects:
		letter = sect[-1]
		sect_len = int(sect[:-1])
		
		if letter == 'M':
			seq2 += seq1[count:count+sect_len] # Add corresponding portion of original seq
			qual2 += qual1[count:count+sect_len]
			count += sect_len

		elif letter == 'D':
			seq2 += '*' * sect_len # Add asterisks for each deleted position relative to the reference sequence
			qual2 += '!' * sect_len # Add ! so that when converted to phred-33 number score will be 0 and won't contribute to consensus score calculation


		elif letter == 'I' or letter == 'S':
			count += sect_len

	return seq2, qual2


def get_read_support(sam_dict, contig, position):
	""" Find the base calls and quality scores from reads that align to a given position in a sam file.
	
		Args:
			sam_dict (dict): dict describing the reads aligned at a given location. Made of nested dicts of format 
							{contig: {location_bin: {read_ID : snpclasses.SAM_data class instance}}}
			contig (str): Name of the contig containing the position of interest
			position (int): Position of interest
	
		Returns:
			tuple: bases (list) base calls of corresponding position in reads aligning at the position of interest,
				   qualities (list) quality scores of corresponding position in reads aligning at the position of interest,
				   ids (list) IDs of reads that align to the position of interest.
	"""
	snp_bin = (position//100)*100

	bases, qualities, ids = [], [], []
	if snp_bin in sam_dict[contig].keys():
		for entry in sam_dict[contig][snp_bin].values():
			for i in entry:
				if i.pos <= position and position < i.end:
					try:
						bases.append(i.seq[position - i.pos])
						qualities.append(ord(i.qual[position - i.pos])-33)
						ids.append(i.qname)
					except Exception as e: 
						print("Error processing read %s in contig %s at position %i: " + str(e) %(i.qname, contig, position - i.pos))

	return bases, qualities, ids


def calc_consensus_score(called_base, read_bases, read_qualities, contig, position):
	"""Calculate a measure of read consensus about the  base call at a goven position in a reference genome
	
		Args:
			called_base (str): the base call given by your variant calling tool
			read_bases (list): The bases at the corresponding location in all reads that align to the position of interest
			read_qualities (list): The quality scores at the corresponding location in all reads that align to the position of interest
			contig (str): Name of contig containing the position of interest
			position (int): Position of interest
	
		Returns:
			str: consensus score = sum of quality scores for bases agreeing with base call / sum of all quality scores for position of interest.
	"""
	agree_score = 0
	disagree_score = 0

	for i in range(len(read_bases)):
		base = read_bases[i]
		quality = int(read_qualities[i])
		if base == called_base:
			agree_score += quality
		else:
			disagree_score += quality

	try:
		consensus_score = agree_score / (agree_score + disagree_score)
	except Exception as e:
		print("Exception while calculating consensus score: %s" % e)
		print(contig)
		print(str(position))
		print(called_base)
		print(read_bases)
		print(read_qualities)
		consensus_score = 0

	return "{:0.2f}".format(consensus_score)


def do_MP_SNP_support(sam_file):
	"""Process the sam file for one of your isolates and write a file summarizing the support for variant calls for this isolate.
	
		Args:
			sam_file (str): Path to sam file of reads of isolate to assess SNP support
	"""
	print("Processing %s..." % sam_file)

	isolate_id = sam_file.split('/')[-1].split('.')[0] # Captures a name for this isolate by taking the file name upto the first . assuming the only . character precedes the extension

	sam_dict = {}
	relevent_contigs = {} # dict to store contigs with snps so we can only import reads aligning to those contigs.

	for snp in snps_list:
		if snp.contig not in relevent_contigs.keys():
			relevent_contigs[snp.contig] = ''

	print("Searching %s for %i snps..." %(sam_file, len(snps_list)))
	
	print("Reading %s sam file..." % sam_file)

	if args.pickle and os.path.isfile(sam_file[:-4] + "_sam_dict.pickle.gzip"): # If you've already pickled the dict of reads mapping to a snp location. 
	# Mostly useful for tweaking settings and rerunning.
		try:
			with gzip.open(sam_file[:-4] + "_sam_dict.pickle.gzip", 'rb') as picklein:
				sam_dict = pickle.load(picklein)
		except Exception as e:
			print("pickle_save_error: " + str(e))

	else:
		with open(sam_file, 'r') as sam:
			for line in sam.readlines():
				if line[0] != "@":
					entry = snpclasses.SAM_data(line)
					try:
						if entry.rname in relevent_contigs.keys():
							if any([x in entry.cigar for x in ['D','I','S']]): # Check if the read need to be padded to account for the alignment to reference.
								entry.mod_seq, entry.mod_qual = cigar_mod_read(entry.seq, entry.qual, entry.cigar)
								entry.refresh()
							span_locs = [x for x in range((entry.pos//100)*100, entry.end, 100)] # identify genome bins read is present in
							if entry.rname in sam_dict.keys(): # Have we already seen this contig?
								for loc in span_locs: # Process each genome bin that this read is present in.
									if loc in sam_dict[entry.rname].keys(): # If we've already seen a read aligned to this bin
										if entry.qname in sam_dict[entry.rname][loc].keys(): # If partner reads have the same name
											if entry.tlen not in [current.tlen for current in sam_dict[entry.rname][loc][entry.qname]]: # Should be positive for one read and negative for the other if they map to the same contig.
												sam_dict[entry.rname][loc][entry.qname].append(entry)
										else:
											sam_dict[entry.rname][loc][entry.qname] = [entry]
									else:
										sam_dict[entry.rname][loc] = {entry.qname:[entry]}
							else:
								sam_dict[entry.rname] = {loc:{entry.qname:[entry]} for loc in span_locs}
					except Exception as e:
						
						print("sam_error: " + str(e))
		if args.pickle:
			try:	
				with gzip.open(sam_file[:-4] + "_sam_dict.pickle.gzip", 'wb') as pickleout:
					pickle.dump(sam_dict, pickleout, protocol=pickle.HIGHEST_PROTOCOL)
			except Exception as e:
				print("pickle_load_error: " + str(e))



	print("Processing %s reads..." % sam_file)

	try:
		print_line_list = []



		for snp in snps_list:
			position = snp.position
			snp_bin = (position//100)*100
			contig = snp.contig


			print_line = [contig, str(position), contig + '_' + str(position)]


			base_calls, base_qualities, read_ids = get_read_support(sam_dict, contig, position)

			try:
				variant_call = max(set(base_calls), key=base_calls.count)
			except:
				print("No reads were found in {} mapping to the snp in contig {} at position {}".format(sam_file, contig, position))
				variant_call = "?"

			ref_base = ref_allele_dict[contig][position]
			
			if args.calls_file:
				if isolate_id in calls_dict.keys():
					if position in calls_dict[isolate_id][contig].keys():
						assembly_base = calls_dict[isolate_id][contig][position]
					else:
						assembly_base = '?'
				else:
					assembly_base = '?'
			else:
				assembly_base = '?'

			if ref_base == variant_call:
				print_line.append('1')
			else:
				print_line.append('0')

			if assembly_base == variant_call: # Base in the assembly might not be the same as the variant call if the variant
			# caller requires a coverage or some other threshold to call a variant.
				print_line.append('1')
			else:
				print_line.append('0')

			print_line += [ref_base, variant_call, assembly_base, str(len(base_calls))]
			
			try:
				print_line.append(str(calc_consensus_score(variant_call, base_calls, base_qualities, contig, position)))
			except Exception as e:
				print("Exception calculating consensus score: " + str(e))
			print_line.append(",".join(base_calls))
			print_line.append(",".join([str(x) for x in base_qualities]))
			print_line.append(",".join(read_ids))

			print_line_list.append(print_line)

	except Exception as e:
		print("Exception processing reads: " + str(e) + " on line " + str(sys.exc_info()[-1].tb_lineno) + " for contig " + str(contig) + ", position " + str(position))

	try:
		if not args.outdir[-1] == '/':
			outdir = args.outdir + '/'
		else:
			outdir = args.outdir

		outfile = outdir + sam_file.split('/')[-1][:-4] + "_snp_support.txt"

		print("Writing %s snp support file to %s..." % (sam_file, outfile))

		with open(outfile, 'w+') as outfile:
			outfile.write("Contig\tPosition\tContig_Position\tReference_Reads_agree\tAssembly_Reads_agree\tReference_base\tIsolate_base\tAssembly_base\tCoverage\tConsensus_score\tRead_base_calls\tRead_qualities\tRead_identifiers\n")
			outfile.write("\n".join(["\t".join(line) for line in print_line_list]) + '\n')
		outfile.close()
	except Exception as e:
		print("Error saving snp_support file: {}".format(e))



def pool_MP_snp_check(sam_files, threads, chunksize):
	"""Manager of worker processes
	
		Args:
			sam_files (list): List of the sam files for isolates you want to process
			threads (int): Number of threads to use
			chunksize (int): approximate number of sam files to give to each process to work through.
	"""
	pool = multiprocessing.Pool(processes=threads)

	pool.imap_unordered(do_MP_SNP_support, sam_files, chunksize)
	pool.close()
	pool.join()



parser = argparse.ArgumentParser(
	description="Given a SNP list file of format contig_ID\tsnp_position, and .sam files, looks up the read support for each SNP in each assembly and reports that information, highlighting any troubling low confidence base-calls.")
parser.add_argument(
	"-snp",  dest="snp_file", required = True,
	help="Specify file containing a list of SNP locations. 1 per line. "
	)
parser.add_argument(
	"-sams",  dest="sam_files", required = True, nargs='+',
	help="Specify all .sam files you want read support for the positions in your snp file. "
	)
parser.add_argument(
	"-ref",  dest="ref_file", required = True,
	help="Specify file with reference genome in fasta fomat (N.B. headers of contigs must match contig names in SNP file). "
	)
parser.add_argument(
	"-calls",  dest="calls_file", required = False,
	help="Optional. Specify file with base calls from your assembler. "
	)
parser.add_argument(
	"-outdir",  dest="outdir", required = True,
	help="Specify output directory. "
	)
parser.add_argument(
	"-threads",  dest="threads", type=int, nargs="?", default = 1,
		help="Specify number of threads to use. Default: 1"
	)
parser.add_argument(
	"-pickle", action='store_true',
		help="Optional storage of reads that map to regions with snps as a pickled dict of format {contig: {1kb_bin: [list of reads]}}"
	)

args = parser.parse_args(sys.argv[1:])

chunksize = len(args.sam_files) // args.threads

snps_list = []
with open(args.snp_file, 'r') as snps_file:
		for snps_line in snps_file.readlines():
			try:
				snps_list.append(snpclasses.snp(snps_line))
			except Exception as e:
				print('Exception processing SNPs in do_MP_SNP_support function: ' + str(e))


with open(args.ref_file, 'r') as ref:
	ref_fasta_dict = fasta_to_dict(ref.read())


ref_allele_dict = {}

for snp in snps_list:
	if snp.contig not in ref_allele_dict.keys():
		ref_allele_dict[snp.contig] = {snp.position : ref_fasta_dict[snp.contig][int(snp.position)-1]} #snp positions are base 1. -1 to convert them to base 0 positions.
	else:
		ref_allele_dict[snp.contig][snp.position] = ref_fasta_dict[snp.contig][int(snp.position)-1]

if args.calls_file:

	calls_dict = {}

	with open(args.calls_file, 'r') as callsf:
		contents = callsf.readlines()
		headers = contents[0].split()[1:]
		for row in contents[1:]:
			cols = row.split()
			calls_dict[cols[0]] = {}
			for i, base in enumerate(cols[1:]):
				contig = headers[i].split('_')[0]
				pos = headers[i].split('_')[1]
				if contig not in calls_dict[cols[0]].keys():
					calls_dict[cols[0]][contig] = {pos : base}
				else:
					calls_dict[cols[0]][contig][pos] = base


pool_MP_snp_check(args.sam_files, args.threads, chunksize)

print("Total run time: %.2f" %(time.time()-start_time))


