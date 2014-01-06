'''
mergeStats.py
for pipe_vda.py

merges the general statistics for two set of reads.

example:
python mergeStats.py <new_run_stats.tab> sdOutgroupMutiplier [reads_to_replace] <mergeDirectory> 

Created:	23/10/2012
Modified:	18/02/2013
author: David Edwards
'''
import sys, glob
from pipe_utils import splitPath

inFileName = sys.argv[1]
(inPrefix, inName, inExt) = splitPath(inFileName)
inFile = open(inFileName)
sdOutgroupMutiplier = int(sys.argv[2])
mergeFileName = sys.argv[4] + inName + inExt
mergeFile = open(mergeFileName, "r")
replace = sys.argv[3]

average = 0
sd = 0
count = 0
output = ""

# Combine the two sets

# First get the merge run list
# setting any reads in the replace list to "fail"
for line in mergeFile:
    items = line.split()
    if items[0] != "Name":
        if items[-1] != "f":
            # if there is a replaceRead list and the read appears on the list
            # change it to a fail...
            if replace != "":
                for name in replace.split(", "):
                    if ("'"+items[0]+"'") == name:
                        line = line[:-2] + "f\n"
                        items = line.split()
            if items[-1] != "f":
                count += 1
                number = float(items[6]) / (float(items[1])/100)
                average += number
                sd += (number*number)
    output += line
mergeFile.close()

# then add the new run stats.tab leaving off the header
for line in inFile:
    items = line.split()
    if items[0] != "Name":
        if items[-1] != "f":
            count += 1
            number = float(items[6]) / (float(items[1])/100)
            average += number
            sd += (number*number)
        output += line
inFile.close()

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
            if number < (average - (sd * sdOutgroupMutiplier)) or number > (average + (sd * sdOutgroupMutiplier)):
                if items[-1] != "o":
                    finalOutput += line[:-1] + "o\n"
                else:
                    finalOutput += line + "\n"
            else:
                if items[-1] != "i":
                    finalOutput += line[:-1] + "i\n"
                else:
                    finalOutput += line + "\n"

mergeFile = open(mergeFileName, "w")
mergeFile.writelines(finalOutput)
mergeFile.close()