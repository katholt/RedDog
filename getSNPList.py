'''
getSNPList.py

takes a number of -q30 vcfs based on the same reference and creates a list of all positions in the reference where a homozygous SNP occurs in at least one of the vcfs.

outputs list to user-defined file

example:
python getSNPList.py <statsFile> <output>

Created:	26/3/2011
Modified:	15/5/2011
author: David Edwards
'''
import sys, glob
from pipe_utils import splitPath

finalList = []
(prefix, name, ext) = splitPath(sys.argv[1])
vcfs = prefix + '/vcf/*_q30.vcf'
for file in glob.glob(vcfs):
    (prefix, name, ext) = splitPath(file)    
    statsFile = open(sys.argv[1])
    for sample in statsFile:
        if sample.startswith("Name") != True:
            splitSample = sample.split("\t")
            if (name[:-4] == splitSample[0]) and (splitSample[-1].startswith("i") == True):
                vcfFile = open(file)
                varList = []
                mergeList = []
                for line in vcfFile:
                    if line.startswith("#") == False:
                        splitLine = line.split("\t")        
                        if splitLine[7].startswith("I") != True:
                            if splitLine[4].find(",") == -1:
                                if splitLine[3] != splitLine[4]:
                                    varList.append(int(splitLine[1]))
                if len(varList) > 0:
                    if len(finalList) < 1:
                        finalList = varList
                    else:
                        posFL = 0
                        posVL = 0
                        lenFL = len(finalList)
                        lenVL = len(varList)
                        while lenFL > 0 and lenVL > 0:
                            if finalList[posFL] <= varList[posVL]:
                                if finalList[posFL] == varList[posVL]:
                                    posVL = posVL + 1
                                    lenVL = lenVL - 1
                                mergeList.append(finalList[posFL])
                                posFL = posFL + 1
                                lenFL = lenFL - 1
                            else:
                                mergeList.append(varList[posVL])
                                posVL = posVL + 1
                                lenVL = lenVL - 1
                        if lenFL == 0 and lenVL > 0:
                            for varPosition in varList[(posVL):]:
                                mergeList.append(varPosition) 
                        elif lenVL == 0 and lenFL > 0:
                            for varPosition in finalList[(posFL):]:
                                mergeList.append(varPosition)
                        finalList = mergeList
        else:
            print "header"
    statsFile.close()
output = ""
for varPosition in finalList:
    output = output + str(varPosition) +"\n"
outputFile = open(sys.argv[2], "w")
outputFile.write(output)
outputFile.close()