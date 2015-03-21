'''
Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)
'''
# input = reference sequence and table of SNP alleles
# output = mfasta whole genome alignment of pseudosequences
import string, re, collections
import os, sys, subprocess
from optparse import OptionParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq

def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-r", "--ref", action="store", dest="ref", help="reference (fasta)", default="")
	parser.add_option("-s", "--snps", action="store", dest="snps", help="snp table (CSV)", default="")
        parser.add_option("-o", "--offset", action="store", dest="offset", help="offset (added to SNP positions in table to match reference coordinates)", default="0")
#	parser.add_option("-f", "--format", action="store", dest="format", help="format (fasta)", default="fasta")

	return parser.parse_args()

if __name__ == "__main__":

	(options, args) = main()
	
	# read in ref seq
#	ref = SeqIO.read(options.ref, options.format)
	ref = SeqIO.read(options.ref, "fasta")

	# read in SNP info
	f = file(options.snps,"r")
	strains = []
	snps = collections.defaultdict(dict) # key 1 = strain, key2 = pos (corrected by offset), value = allele
	for line in f:
		fields = line.rstrip().split(",")
		if len(strains) > 0:
			pos = int(fields[0]) + int(options.offset)
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
		print ">" + strain
		print str(seq)
		
