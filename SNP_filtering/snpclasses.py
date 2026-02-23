#!/usr/bin/env python3

# AUTHOR      :  ALAN COLLINS
# VERSION     :  v1
# DATE        :  2021-1-8
# DESCRIPTION :	 classes for processing sam and snp location files.

class SAM_data(object):
	"""stores columns of SAM entry as attributes
	SAM format columns:

Col		Field		Type		Regexp/Range					Brief description
1 		QNAME 		String 		[!-?A-~]{1,254} 				Query template NAME
2 		FLAG 		Int 		[0, 2^16 − 1] 					bitwise FLAG
3 		RNAME 		String 		\*|[:rname:∧*=][:rname:]* 		Reference sequence NAME11
4 		POS 		Int 		[0, 2^31 − 1] 					1-based leftmost mapping POSition
5 		MAPQ 		Int 		[0, 2^8 − 1] 					MAPping Quality
6 		CIGAR 		String 		\*|([0-9]+[MIDNSHPX=])+ 		CIGAR string
7 		RNEXT 		String 		\*|=|[:rname:∧*=][:rname:]*		Reference name of the mate/next read
8 		PNEXT 		Int 		[0, 2^31 − 1] 					Position of the mate/next read
9 		TLEN 		Int 		[−2^31 + 1, 2^31 − 1] 			observed Template LENgth
10 		SEQ 		String 		\*|[A-Za-z=.]+ 					segment SEQuence
11 		QUAL 		String 		[!-~]+ 							ASCII of Phred-scaled base QUALity+33

original read len = sum of numbers in cigar string. e.g. 56H147M means read is 56+147 long = 203.
"""
	def __init__(self, object):
		self.qname = object.split('\t')[0]
		self.flag = object.split('\t')[1]
		self.rname = object.split('\t')[2]
		self.pos = int(object.split('\t')[3])
		self.mapq = int(object.split('\t')[4])
		self.cigar = object.split('\t')[5]
		self.rnext = object.split('\t')[6]
		self.pnext = object.split('\t')[7]
		self.tlen = object.split('\t')[8]
		self.seq = object.split('\t')[9]
		self.qual = object.split('\t')[10]
		self.ln = len(self.seq)
		self.end = self.pos + self.ln
		self.mod_seq = ''
		self.mod_qual = ''
		self.refreshed = False

	def refresh(self):
		if not self.refreshed:
			self.seq = self.mod_seq
			self.qual = self.mod_qual
			self.ln = len(self.seq)
			self.end = self.pos + self.ln
			self.refreshed = True

class snp():
	"""Simple class to store contig ID and position info for SNPs"""
	def __init__(self, line):
		self.contig = line.split()[0].split(':')[0]
		self.position = int(line.split()[1])
		


