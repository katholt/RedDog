'''
collateRepStats.py
for mapping_pipe.py


collates the general statistics for set of reads that have gone through
the pipeline during a pangenome or phylogeny run for a particular replicon

example:
python collateRepStats.py reference_name <example>_rep_cover.txt replicon sd_multiplier runType

Created:	29042013
Modified:	28102013 changed to mean + SD * sd_multiplier only
author: David Edwards
'''
import sys, glob
from pipe_utils import (splitPath)

reference_name = sys.argv[1]
example_rep_cover_File = sys.argv[2]
replicon = sys.argv[3]
sd_multiplier = float(sys.argv[4])
runType = sys.argv[5]

(prefix, middle, ext) = splitPath(example_rep_cover_File)
outfile_RepStats_name = prefix[:-4] + reference_name + '_' + replicon +  '_RepStats.txt'

output_RepStats = ''
header_RepStats = 'Isolate\tCover%_' + replicon + '\tDepth_' + replicon + '\tMapped%_' + replicon + '\tMapped%_Total\tTotal_Reads\tSNPs\tHets_Removed\tIndels\tIngroup/Fail\n'
inFiles = prefix + '/*_' + replicon +'_RepStats.txt'
if runType == 'pangenome':
	for inFile in glob.glob(inFiles):
	    repStatsFile = open(inFile)
	    output_RepStats += repStatsFile.read()
else:
	average = 0
	sd = 0
	count = 0
	for inFile in glob.glob(inFiles):
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
	for file in glob.glob(inFiles):
	    statsFile = open(file)
	    for line in statsFile:
	        stats = line.split()
	        output_RepStats += line[:-1]
	        if stats[-1] != "f":
	            number = float(stats[6]) / (float(stats[1])/100)
#	            if number < (average - (sd * sd_multiplier)) or number > (average + (sd * sd_multiplier)):
	            if number > (average + (sd * sd_multiplier)):
	                output_RepStats += "\to\n"
	            else:
	                output_RepStats += "\ti\n"
	        else:
	            output_RepStats += "\n"
	    statsFile.close()

outfile_RepStats = open(outfile_RepStats_name,"w")
outfile_RepStats.write(header_RepStats)
outfile_RepStats.writelines(output_RepStats)
outfile_RepStats.close()
