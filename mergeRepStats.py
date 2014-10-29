'''
mergeRepStats.py
for RedDog.py

merges the general 'replicon' statistics for two set of reads.

example:
python mergeRepStats.py <new_replicon_RepStats.tab> sdOutgroupMutiplier reads_to_replace <mergeDirectory> runType 

Created:	23/10/2012
Modified:	28/10/2013 to mergeRepStats from mergeStats
            15/04/2014 changed to produce outgroup.txt file if there are any outgroups to report
            20/05/2014 fix to outgroup reporting
            04/07/2014 change to replace_reads handling

author: David Edwards
'''
import sys, glob
from pipe_utils import splitPath

inFileName = sys.argv[1]
(inPrefix, inName, inExt) = splitPath(inFileName)
sdOutgroupMutiplier = int(sys.argv[2])
mergeFileName = sys.argv[4] + inName + inExt
replace = sys.argv[3]
runType = sys.argv[5]
outgroup_outfile_name = sys.argv[4] + inName[:-9] +  '_outgroups.txt'
outgroups = []

average = 0
sd = 0
count = 0
output = ""

# Combine the two sets
# note: only need to count SNPs for 'phylogeny' runType

# First get the merge run list
# setting any reads in the replace list to "fail"
mergeFile = open(mergeFileName, "r")
for line in mergeFile:
    items = line.split()
    if items[0] != "Isolate":
        if items[-1] != "f":
            # if there is a replaceRead list and the read appears on the list
            # change it to a fail...
            if replace != "":
                for name in replace.split(","):
                    if items[0] == name:
                        line = line[:-2] + "f\n"
                        items = line.split()
            if items[-1] != "f" and runType=='phylogeny':
                count += 1
                number = float(items[6]) / (float(items[1])/100)
                average += number
                sd += (number*number)
    output += line
mergeFile.close()

# then add the new run stats.tab leaving off the header
inFile = open(inFileName)
for line in inFile:
    items = line.split()
    if items[0] != "Isolate":
        if items[-1] != "f" and runType=='phylogeny':
            count += 1
            number = float(items[6]) / (float(items[1])/100)
            average += number
            sd += (number*number)
        output += line
inFile.close()

# work out new mean and sd
# Note: 'pangenome' runType will 'skip' this as count == 0 
if count > 0:
    average = average/(count*1.0)
    sd = ((sd*1.0)/count - average*average)
    if sd < 0:
        sd = sd * -1
    sd = sd**(1/2.0)

# identify ingroups/outgroups in the merged set for 'phylogeny' runType
finalOutput = ""
for line in output.split("\n"):
    if line != "":
        items = line.split()
        if items[0] == "Isolate" or items[-1] == "f":
            finalOutput += line + "\n"
        elif runType == 'phylogeny':
            number = float(items[6]) / (float(items[1])/100)
#            if number < (average - (sd * sdOutgroupMutiplier)) or number > (average + (sd * sdOutgroupMutiplier)):
            if number > (average + (sd * sdOutgroupMutiplier)):
                outgroups.append(items[0])
                if items[-1] != "o":
                    finalOutput += line[:-1] + "o\n"
                else:
                    finalOutput += line + "\n"
            else:
                if items[-1] != "i":
                    finalOutput += line[:-1] + "i\n"
                else:
                    finalOutput += line + "\n"
        else:
            if items[-1] != "i":
                finalOutput += line[:-1] + "i\n"
            else:
                finalOutput += line + "\n"

if outgroups != []:
    outgroup_outfile = open(outgroup_outfile_name,"w")
    for outgroup in outgroups:
        outgroup_outfile.write(outgroup+'\n')
    outgroup_outfile.close()

mergeFile = open(mergeFileName, "w")
mergeFile.writelines(finalOutput)
mergeFile.close()