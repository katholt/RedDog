'''
collateStats.py
for mapping_pipe.py


collates the general statistics for set of reads that have gone through
the pipeline

example:
python collateStats.py <inputPath> minDepth sdOutgroupMutiplier <ref_tmap_stats.tab>

Created:	2/3/2012
Modified:	18/02/2013
author: David Edwards
'''
import sys, glob
outFile = open(sys.argv[4], "w")
inFiles = sys.argv[1] + "*_stats.txt"
output = ""
average = 0
sd = 0
count = 0
header = "Name" +"\t"
minDepth = sys.argv[2]
for file in glob.glob(inFiles):
    statsFile = open(file)
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
sdOutgroupMutiplier = int(sys.argv[3])
for file in glob.glob(inFiles):
    statsFile = open(file)
    for line in statsFile:
        stats = line.split()
        output += line[:-1]
        if stats[-1] != "f":
            number = float(stats[6]) / (float(stats[1])/100)
            if number < (average - (sd * sdOutgroupMutiplier)) or number > (average + (sd * sdOutgroupMutiplier)):
                output += "\to\n"
            else:
                output += "\ti\n"
        else:
            output += "\n"
    statsFile.close()
header = header + "%_Reference\tDepth\tDepth>=" + minDepth +"\tTotal_Reads\t%_Mapped\tq30_SNPs\tq30_hets_removed\tq30_INDELs\tin/outgroup/fail\n"
outFile.write(header)
outFile.writelines(output)
outFile.close()