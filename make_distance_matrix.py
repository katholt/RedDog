#!/usr/bin/env python
#
# make_distance_matrix.py
#
# - options version
#
# Calculates the distance matrix based on the pair-wise differences in SNPS
# in the allele matrix
#
# 
# Outputs the file to NEXUS format
#
# Example command:
'''
module load python-gcc/2.7.5
python make_distance_matrix.py -i my_strains_alleles.csv
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
	parser.add_option("-o", "--output_file", action="store", dest="output_file", help="name for output file (no ext, only required if '_alleles' missing from input name), Note: if '_alleles(_*)' ending present, this is removed and '_SNP_diff.nxs' added [setting for RedDog pipeline] (default: none)", default="")
	parser.add_option("-d", "--directory", action="store", dest="directory", help="directory to send output files (default: none)", default="")
		
	return parser.parse_args()

if __name__ == "__main__":

	(options, args) = main()

	### MAIN PROCESS
	inputName = options.input_file
	inputFile = open(inputName)

	(prefix, name, ext) = splitPath(inputName)

	strains = []

	if name.find('_alleles') != -1:
		split_name = name.split('_')
		split = 0
		refName = split_name[split]
		split += 1
		while split_name[split] != 'alleles':
			refName += '_' + split_name[split]
			split += 1 
		directory = options.directory		
		if directory == '':
			if prefix != '':
				outName = prefix + '/' + refName + '_SNP_diff.nxs'
			else:
				outName = refName + '_SNP_diff.nxs'
		elif directory.endswith('/'):
			outName = directory + refName + '_SNP_diff.nxs'
		else:
			outName = directory + '/' + refName + '_SNP_diff.nxs'
	else:
		if options.output_file != '':
			strains.append(options.reference)
			directory = options.directory
			if directory == '':
				if prefix != '':
					outName = prefix + '/' + name + '_SNP_diff.nxs'
				else:
					outName = name + '_SNP_diff.nxs'
			elif directory.endswith('/'):
				outName = directory + name + '_SNP_diff.nxs'
			else:
				outName = directory + '/' + name + '_SNP_diff.nxs'
		else:
			print "No output file specified: check option -o output_file"
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
	for count1 in range(1,len(strains)):
		differenceCount.append([0])
		while len(differenceCount[count1-1]) < count1:
			differenceCount[count1-1].append(0)

	line = inputFile.readline()
	while line != "":
		bases = line.rstrip().split(',')
		for count1 in range(2,(len(bases)-1)):
			for count2 in range((count1+1),len(bases)):
				if bases[count1] != bases[count2] and bases[count1] != '-' and bases[count2] != '-':
					differenceCount[count2-3][count1-2] += 1
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
			output = output + str(differenceCount[count1-1][count2-1])
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
