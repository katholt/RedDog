'''
collateMergeStats.py

Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)

merges the general statistics for merged set of reads.

example:
python collateMergeStats.py <run_stats.tab> sdOutgroupMutiplier "mergedReadsToFail" <mergeStatsTxtDirectory> 

Created:	02/11/2012
Modified:	18/02/2013
author: David Edwards
'''
import sys, glob
from pipe_utils import splitPath

fileName = sys.argv[1]
sdOutgroupMultiplier = int(sys.argv[2])
statsFilesPattern = sys.argv[4] + "*_merged_stats.txt"
statsFiles = []
if type(statsFilesPattern) == list:
    for pattern in statsFilesPattern:
        statsFiles.append(glob.glob(pattern))
else:
    statsFiles = glob.glob(statsFilesPattern)

replace = sys.argv[3]

average = 0
sd = 0
count = 0
output = ""

# Combine the two sets

# First get the run stats.tab
# setting any reads in the replace list to "fail"
inFile = open(fileName, "r")
for line in inFile:
    items = line.split()
    if items[0] != "Name":
        if items[-1] != "f":
                # if the read appears on the list
                # change it to a fail...
            if replace != "-":
                for name in replace.split(" "):
                    if (items[0]) == name:
                        line = line[:-2] + "f\n"
                        items = line.split()
            if items[-1] != "f":
                count += 1
                number = float(items[6]) / (float(items[1])/100)
                average += number
                sd += (number*number)
    output += line
inFile.close()


# then add the merged bams stats

for files in range(len(statsFiles)):
    statsFile = open(statsFiles[files], "r")
    for line in statsFile:
        items = line.split()
        if items[-1] != "f":
            count += 1
            number = float(items[6]) / (float(items[1])/100)
            average += number
            sd += (number*number)
        output += line
    statsFile.close()

# work out new mean and sd
if count > 0:
    average = average/(count*1.0)
    sd = ((sd*1.0)/count - average*average)
    if sd < 0:
        sd = sd * -1
    sd = sd**(1/2.0)

# identify ingroups/outgroups in the merged set
finalOutput = ""
for line in output.split("\n"):
    if line != "":
        items = line.split()
        if items[0] == "Name" or items[-1] == "f":
            finalOutput += line + "\n"
        else:
            number = float(items[6]) / (float(items[1])/100)
            if number < (average - (sd * sdOutgroupMultiplier)) or number > (average + (sd * sdOutgroupMultiplier)):
                if items[-1] != "o":
                    if items[-1] == "i":
                        finalOutput += line[:-1] + "o\n"
                    else:
                        finalOutput += line + "\to\n"    
                else:
                    finalOutput += line + "\n"
            else:
                if items[-1] != "i":
                    if items[-1] == "o":                    
                        finalOutput += line[:-1] + "i\n"
                    else:
                        finalOutput += line + "\ti\n"
                else:
                    finalOutput += line + "\n"

outFile = open(fileName, "w")
outFile.writelines(finalOutput)
outFile.close()