'''
make_distance_matrix.py

- version for use with pipeline

Calculates the distance matrix based on the pair-wise differences in SNPS
in the allele matrix

Outputs the file to NEXUS format

Examples:
python make_distance_matrix.py <allele matrix>

python make_distance_matrix.py *_alleles.csv

Created:	10/02/2013
Modified:	25/02/2013
author: David Edwards
'''
import sys, glob
from pipe_utils import splitPath
inputName = sys.argv[1]
inputFile = open(inputName)
(prefix, name, ext) = splitPath(inputName)
refName = name[:-8]
strains = []
strains.append(refName)
outName = prefix + '/' + refName + '_SNP_diff.nxs'

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
					differenceCount[count1-1][count2-count1-1] += 1
	line = inputFile.readline()

longestName = 0
for strain in range(len(strains)):
	if longestName <= len(strains[strain]):
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
