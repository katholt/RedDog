#!/usr/bin/env python3

'''
Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)

Updated to python3 by Kelly Wyres, 2018
'''



# input = reference sequence and table of SNP alleles
# output = mfasta whole genome alignment of pseudosequences
import string, re, collections
import os, sys, subprocess
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq

def get_arguments():
	parser = argparse.ArgumentParser(description='create whole genome pseudo aln using snp matrix and reference')
	parser.add_argument("-r", "--ref", action="store", dest="ref", help="reference (fasta)", default="")
	parser.add_argument("-s", "--snps", action="store", dest="snps", help="snp table (CSV)", default="")
	parser.add_argument("-o", "--offset", action="store", dest="offset", help="offset (added to SNP positions in table to match reference coordinates)", default="0")
	
	return parser.parse_args()

def main():

	# Get arguments
	args = get_arguments()
	
	# read in ref seq
	ref = SeqIO.read(args.ref, "fasta")

	# read in SNP info
	f = open(args.snps,"r")
	strains = []
	snps = collections.defaultdict(dict) # key 1 = strain, key2 = pos (corrected by offset), value = allele
	for line in f:
		fields = line.rstrip().split(",")
		if len(strains) > 0:
			pos = int(fields[0]) + int(args.offset)
			for i in range(1,len(fields)):
				snps[strains[i]][pos] = fields[i]
		else:
			strains = fields
			
	# print out mutated copy of sequence for each strain
	for strain in snps:
		# make mutable copy of reference seq
		seq = ref.seq.tomutable()
		for snp in snps[strain]:
			seq[snp-1] = snps[strain][snp]
		print(">" + strain)
		print(str(seq))
		
if __name__ == '__main__':
	main()