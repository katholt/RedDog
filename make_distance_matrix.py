#!/usr/bin/env python
#
# make_distance_matrix.py
#
# - options version
#
# Calculates the distance matrix based on the pair-wise differences in SNPS
# in the allele matrix
#
# default settings will produce matrix based on the pair-wise differences in SNPS
# each difference (SNP) increases pair-wise difference score by 1
#
# Extended difference models
# (Note: only use these if all SNPs have only one variant - script cannot
# handle situation where more than one variant occurs at a SNP position)
# If any of the following three are activated, the pair-wise difference score 
# for the SNP (1) is multiplied by the following consequence(s), which are
# cumulative in effect if more than one test specified: 
#
# -t includes test for whether SNP is transition or transversion
#		consequence: transition - 0.5, transversion - 1
# -c includes test for whether SNP is in coding or noncoding region of genome
#		consequence: noncoding - 0.5, coding - 1
# -s includes test for whether coding SNP is synonymous or nonsynonymous 
#		consequence: synonymous - 1, nonsynonymous - 2
# 
# For example, if -t -c is specified and a SNP is a coding transversion, the
# pair-wise difference score will be 1 (1*1*1), whilst a noncoding transition
# will score 0.25 (1*0.5*0.5)
#
# Using either -c or -s requires the SNP consequences file from parseSNPtable
#
# Outputs the file to NEXUS format
#
# Example command:
'''
module load python-gcc/2.7.5
python make_distance_matrix.py -i my_strains_alleles.csv -o my_strains_SNP_diff.nxs -t -c -s -C my_strains_consequences.txt
'''
#author: David Edwards
#
#Created:	10/02/2013
#Changes:	20/04/2014 -optparse added, extended difference models added

import sys, glob
from pipe_utils import splitPath
from optparse import OptionParser

def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i", "--input_file", action="store", dest="input_file", help="allele table (required)", default="")
	parser.add_option("-o", "--output_file", action="store", dest="output_file", help="name for output file (no ext), otherwise removes '_alleles(_*)' ending and adds '_SNP_diff.nxs' [setting for RedDog pipeline] (default: none)", default="")
	parser.add_option("-r", "--reference", action="store", dest="reference", help="name of reference (default: none, otherwise removes '_alleles(_*)' ending from input file and uses this as reference name [setting for RedDog pipeline] (default: none)", default="")
	parser.add_option("-d", "--directory", action="store", dest="directory", help="directory to send output files (default: none)", default="")

	parser.add_option("-t", "--transition_test", action="store_true", dest="transition_test", help="include transition/transversion test in difference calulation (default: off [False])", default=False)
	
	parser.add_option("-C", "--consequences_file", action="store", dest="consequences_file", help="SNP consequences file: required for coding test (-c) and/or synonymous testing (-s) (default: none)", default="")

	parser.add_option("-c", "--coding_test", action="store_true", dest="coding_test", help="include coding/noncoding test in difference calulation (default: off [False])", default=False)
	parser.add_option("-s", "--synonymous_test", action="store_true", dest="synonymous_test", help="include synonymous/nonsynonymous test (for coding SNPs) in difference calulation (default: off [False])", default=False)
		
	return parser.parse_args()

if __name__ == "__main__":

	(options, args) = main()

	def getTransitionScore(ref_call,alt_call):
		if ((ref_call == 'T' and alt_call == 'C') or (ref_call == 'C' and alt_call == 'T') 
			or (ref_call == 'A' and alt_call == 'G') or (ref_call == 'G' and alt_call == 'A')):
			return 0.5;
		else:
			return 1;

	### MAIN PROCESS
	inputName = options.input_file
	inputFile = open(inputName)

	(prefix, name, ext) = splitPath(inputName)

	strains = []

	if names.find('_alleles') != -1:
		split_name = name.split('_')
		split = 0
		refName = split_name[split]
		split += 1
		while split_name[split] != 'alleles':
			refName += '_' + split_name[split]
			split += 1 
		strains.append(refName)
		directory = options.directory		
		if directory == '':
			outName = prefix + '/' + refName + '_SNP_diff.nxs'
		elif directory.endswith('/'):
			outName = directory + refName + '_SNP_diff.nxs'
		else:
			outName = directory + '/' + refName + '_SNP_diff.nxs'
	else:
		if options.output_file != '' and options.reference != '':
			strains.append(options.reference)
			directory = options.directory
			if directory == '':
				outName = prefix + '/' + refName + '_SNP_diff.nxs'
			elif directory.endswith('/'):
				outName = directory + refName + '_SNP_diff.nxs'
			else:
				outName = directory + '/' + refName + '_SNP_diff.nxs'
		elif options.output_file == '' and options.reference != '':
			print "No output file specified: check option -o output_file"
			sys.exit()
		elif options.output_file != '' and options.reference == '':
			print "No reference specified: check option -r reference"
			sys.exit()
		elif options.output_file == '' and options.reference == '':
			print "No output file specified: check option -o output_file"
			print "No reference specified: check option -r reference"
			sys.exit()

	positions = []
	changes = []
	if coding_test or synonymous_test:
		if options.consequences_file == "":
			print "No consequences file specified: check option -C consequences_file"
			sys.exit()
		else:
			try:
				consequences_file = open(options.consequences_file)
				for line in consequences_file:
					if not line.startswith('SNP') and line != '':
						data = line.split('\t')
						positions.append(data[0])
						changes.append(data[3][0])
			except:
				print "Could not open consequences file: check option -C consequences_file"
				sys.exit()

	line = inputFile.readline()
	labels = line.split(',')
	for label in range(0,len(labels)):
	    if labels[label] != 'Pos' and labels[label] != 'Ref':
	    	if labels[label].find('\n') == -1:
	    		strains.append(labels[label])
	    	else:
	    		strains.append(labels[label][:-1])

	differenceCount = []
	for count1 in range(0,(len(strains) - 1)):
		differenceCount.append([0])
		while len(differenceCount[count1]) < (len(strains) - count1 - 1):
			differenceCount[count1].append(0)

	line = inputFile.readline()
	while line != "":
		bases = line.split(',')
		for count1 in range(1,(len(bases)-1)):
			for count2 in range((count1+1),len(bases)):
				if bases[count1][0:1] != bases[count2][0:1] and bases[count1][0:1] != '-' and bases[count2][0:1] != '-':
					if not transition_test and not coding_test and not synonymous_test:
						differenceCount[count1-1][count2-count1-1] += 1
					else:
						difference = 1
						if transition_test:
							difference *= getTransitionScore(bases[count1][0:1],bases[count2][0:1])
						if coding_test:
							if changes[positions.index(bases[0])] == 'i':
								difference *= 0.5
						if synonymous_test:
							if changes[positions.index(bases[0])] == 'n':
								difference *= 2
						differenceCount[count1-1][count2-count1-1] += difference
		line = inputFile.readline()

	longestName = 0
	for strain in range(len(strains)):
		if longestName < len(strains[strain]):
			longestName = len(strains[strain])


	output = "#NEXUS\n\nBEGIN taxa;\n    DIMENSIONS ntax=" + str(len(strains)) + ";\n"

	output += "TAXLABELS\n    " 
	for strain in range(len(strains)):
		name = strains[strain]
		output = output + strains[strain] + "\n"
		if (strain + 1) < len(strains):
			output += "    "
		else:
			output += ";\nEND;\n\n"

	output = output + "BEGIN distances;\n    DIMENSIONS ntax=" + str(len(strains)) + ";\n"
	output = output + "    FORMAT\n        TRIANGLE=LOWER\n        NODIAGONAL\n    ;\n"
	output = output + "MATRIX\n"

	name = strains[0]
	while len(name) < longestName:
		name += " "
	output = output + name + " \n"
	for count1 in range(1,len(strains)):
		name = strains[count1]
		while len(name) < (longestName + 1):
			name += " "
		output = output + name
		for count2 in range(1,count1+1):
			output = output + str(differenceCount[count2-1][count1-count2])
			if count2 < count1:
				output += "\t"
			else:
				if (count1 + 1) < len(strains):
					output += "\n"
				else:
					output += ";\nEND;\n"

	inputFile.close()
	outputFile = open(outName, "w")
	outputFile.write(output)
	outputFile.close()
