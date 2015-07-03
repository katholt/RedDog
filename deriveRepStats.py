#!/bin/env python
'''
dervieRepStats.py

Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)

collates the general statistics by replicon for set of reads that have gone through
the pipeline - pangenome run (largest replicon, or user-defined list)

example:
python deriveRepStats.py <isolate>_rep_cover.txt replicon depth_fail cover_fail runType mapped_fail check_reads_mapped

Created:	29042013
Modified:	16072014
author: David Edwards
'''
import sys
from pipe_utils import splitPath

output_RepStats = ""

repCoverFile = sys.argv[1]
repCover = open(repCoverFile, "r") 
(prefix, middle, ext) = splitPath(repCoverFile)
seq_name = middle[:-10]

output_RepStats += seq_name +"\t"

replicon = sys.argv[2]
if replicon.find('.') != -1:
    temp_rep = replicon.split('.')
    replicon = temp_rep[0]
replicon_names = []

#get % cover and depth for each replicon
depth_test_value = 0.0
for line in repCover:
	entry = line.split()
	replicon_names.append(entry[0])
	if entry[0] == replicon:
		output_RepStats += entry[3] +"\t"+ entry[2] +"\t"
		depth_test_value = float(entry[2])
		cover_test_value = float(entry[3])

#get % mapped stats for each replicon and the total
replicons_mapped = []
name = prefix + "/" + seq_name + "_samStats.txt"
repMapped = open(name)
for line in repMapped:
	if line.startswith('reads'):
		entry = line.split()
		total_reads = entry[1]
	elif line.startswith('mapped reads'):
		entry = line.split()
		mapped_reads = entry[2]
	elif line.startswith('%'):
		for rep in replicon_names:
			if line.startswith('%'+rep):
				entry = line.split()
				replicons_mapped.append([rep, entry[1]])

mapped_rep = 0.0
for rep in replicon_names:
	value = ""
	for mapped in replicons_mapped:
		if rep == mapped[0]:
			value = mapped[1]
	if rep == replicon:
		if value == "":
			output_RepStats += '0.0\t'			
		else:
			mapped_rep = float(value)*int(mapped_reads)/int(total_reads)
			output_RepStats += str(mapped_rep) +"\t"

mapped_test_value = int(mapped_reads)*100.0/int(total_reads)
output_RepStats += str(mapped_test_value) +"\t"+ total_reads + "\t"

name = prefix + "/q30VarFilter/hets/" + seq_name + "_" + replicon + "_het.txt"
hetFile = open(name)
het = hetFile.read()
hetFile.close()

name = prefix + "/getVCFStats/" + seq_name + "_" + replicon + "_vcf.txt"
vcfFile = open(name)
vcfList = vcfFile.readline()
vcfFile.close()
vcf = vcfList.split()

output_RepStats += vcf[0] +"\t"+ het +"\t"+ vcf[2]

# decide if a sample (strain) is a fail or pass

depthFail = int(sys.argv[3])
coverFail = int(sys.argv[4])
runType = sys.argv[5]
mappedFail = int(sys.argv[6])
check_reads_mapped = sys.argv[7]

line_end = False
if (depth_test_value < depthFail):
	output_RepStats += "\tf\n"
	line_end = True
elif (cover_test_value < coverFail):
	output_RepStats += "\tf\n"
	line_end = True
elif (check_reads_mapped != "off"):
	if check_reads_mapped == replicon:
		if (mapped_test_value < int(mappedFail)):
			output_RepStats += "\tf\n"
			line_end = True
	elif (check_reads_mapped.find(replicon) != -1):
		found_x = False
		final_ratio = 1.0
		list_of_replicons = []
		ratio_of_replicons = []
		check_reps = check_reads_mapped.split(',')
		for item in check_reps:
			if item != 'x':
				if found_x != True:
					list_of_replicons.append(item)
				else:
					ratio_of_replicons.append(float(item))
					final_ratio -= float(item)
			else:
				found_x = True
		ratio_of_replicons.append(final_ratio)
		for i in range(0, len(list_of_replicons)):
			if list_of_replicons[i] == replicon:
				if (mapped_test_value < int(mappedFail)*ratio_of_replicons[i]):
					output_RepStats += "\tf\n"
					line_end = True

if line_end == False:
	if runType == "phylogeny":
		output_RepStats += "\n"
	else:
		output_RepStats += "\ti\n"

name = prefix + "/deriveRepStats/" + seq_name + "_" + replicon + "_RepStats.txt"
out_RepStats = open(name, "w")
out_RepStats.write(output_RepStats)
out_RepStats.close()
