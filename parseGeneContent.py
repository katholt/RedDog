'''
Copyright (c) 2015, David Edwards, Kat Holt
All rights reserved. (see README.txt for more details)
'''
# read in coverage output file (% coverage for each gene in mapped pan genome)
# remove genes that are not covered (to set % level, default 95%) in any strains
# generate presence/absence table (1/0, based on set % level, default 95%)
# generate summary of genes, reporting total number of strains with the gene
# FUTURE:
# could also provide strain IDs and groups, and summarize presence by group
# could also calculate pan and core rarefaction curve data
# could also accept depth file and use this along with % coverage to call presence/absence
# could also accept RAST annotation file and output product identifiers etc along with genes (esp in gene summary table)
#
# Authors: Kat Holt, David Edwards
# created: ?
# last modified: 20141020
#
import string, re
import os, sys, subprocess
from optparse import OptionParser

def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-g", "--gcfile", action="store", dest="gcfile", help="gene content file (csv)", default="")
	parser.add_option("-c", "--cutoff", action="store", dest="cutoff", help="cutoff coverage level (%) (95)", default="95")
	parser.add_option("-o", "--out", action="store", dest="out", help="file to output presence/absence matrix (csv; not output if not provided)", default="")
	parser.add_option("-s", "--summary", action="store", dest="summary", help="output for genewise summary (csv; not output if not provided)", default="")
	
	return parser.parse_args()


if __name__ == "__main__":

	(options, args) = main()
			
	header = []
	if options.out !="":
		o = file(options.out,"w")
	if options.summary !="":
		s = file(options.summary,"w")
		s.write("gene,strainCount\n")
	if options.summary =="" and options.out=="":
		DoError("No output requested. Use -o to generate presence/absence matrix and/or -s to generate gene summary file.")
		
	f = file(options.gcfile,"r")
	for line in f:
		if len(header)==0:
			header = line.rstrip().split(",")
			o.write(line)
		else:
			fields = line.rstrip().split(",")
			gene = fields[0]
			cov = []
			count = 0
			for i in range(1,len(fields)):
				cov.append(float(fields[i]))
				if (float(fields[i]) >= float(options.cutoff)):
					count += 1
			if (count > 0):
				# otherwise exclude gene from outputs
				if options.out !="":
					# convert to presence/absence
					output = fields[0]
					for i in range(1,len(fields)):
						if (float(fields[i]) >= float(options.cutoff)):
							# gene is present
							output += ',1'
						else:
							# gene is absent
							output += ',0'
					output += '\n'
					o.write(output)
				if options.summary !="":
					s.write(gene + "," + str(count) + "\n")
	f.close()
	o.close()
	s.close()
				
				