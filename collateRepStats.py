'''
collateRepStats.py

Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)

collates the general statistics for set of reads that have gone through
the pipeline during a pangenome or phylogeny run for a particular replicon

example:
python collateRepStats.py reference_name prefix replicon sd_multiplier runType sequences_string

Created:	29042013
Modified:	16072014
author: David Edwards
'''
import sys, glob
from pipe_utils import (splitPath)

reference_name = sys.argv[1]
prefix = sys.argv[2]
replicon = sys.argv[3]
sd_multiplier = float(sys.argv[4])
runType = sys.argv[5]
sequences_string = sys.argv[6]
sequences = sequences_string.split(',')

outfile_RepStats_name = prefix + reference_name + '_' + replicon +  '_RepStats.txt'
outgroup_outfile_name = prefix + reference_name + '_' + replicon +  '_outgroups.txt'
outgroups = []

output_RepStats = ''
header_RepStats = 'Isolate\tCover%_' + replicon + '\tDepth_' + replicon + '\tMapped%_' + replicon + '\tMapped%_Total\tTotal_Reads\tSNPs\tHets_Removed\tIndels\tIngroup/Fail\n'
inFiles = []
for sequence in sequences:
    if sequence != '':
        inFiles.append((prefix + 'temp/' + sequence + '/deriveRepStats/' + sequence + '_' + replicon + '_RepStats.txt'))

if runType == 'pangenome':
	for inFile in inFiles:
	    statsFile = open(inFile)
	    output_RepStats += statsFile.read()
	    statsFile.close()
else:
	average = 0
	sd = 0
	count = 0
	for inFile in inFiles:
	    statsFile = open(inFile)
	    for line in statsFile:
	        stats = line.split()
	        if stats[-1] != "f":
	            count += 1
	            number = float(stats[6]) / (float(stats[1])/100)
	            average += number
	            sd += (number*number)
	    statsFile.close()
	if count > 0:
	    average = average/(count*1.0)
	    sd = ((sd*1.0)/count - average*average)
	    if sd < 0:
	        sd = sd * -1
	    sd = sd**(1/2.0)
	for inFile in inFiles:
	    statsFile = open(inFile)
	    for line in statsFile:
	        stats = line.split()
	        output_RepStats += line[:-1]
	        if stats[-1] != "f":
	            number = float(stats[6]) / (float(stats[1])/100)
#	            if number < (average - (sd * sd_multiplier)) or number > (average + (sd * sd_multiplier)):
	            if number > (average + (sd * sd_multiplier)):
	                output_RepStats += "\to\n"
	                outgroups.append(stats[0])
	            else:
	                output_RepStats += "\ti\n"
	        else:
	            output_RepStats += "\n"
	    statsFile.close()

if outgroups != []:
	outgroup_outfile = open(outgroup_outfile_name,"w")
	for outgroup in outgroups:
		outgroup_outfile.write(outgroup+'\n')
	outgroup_outfile.close()

outfile_RepStats = open(outfile_RepStats_name,"w")
outfile_RepStats.write(header_RepStats)
outfile_RepStats.writelines(output_RepStats)
outfile_RepStats.close()
