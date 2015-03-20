'''
collateAllRepGeneCover.py

Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)

Takes a coverage and average depth output divided by replicon__gene, and named 
after each isolate, and converts them to a coverage matrix and average depth matrix,
divided by replicon__gene and strain.

example:
python collateAllGeneCover.py <inputDir> <outputDir> refName sequences_string

Created:	27/5/2013
Modified:   21/07/2014
author: David Edwards
'''
import sys, glob
from pipe_utils import splitPath

inPath = sys.argv[1]
outPath = sys.argv[2]
refName = sys.argv[3]
sequences_string = sys.argv[4]
sequences = sequences_string.split(',')

coverFiles = []
for sequence in sequences:
    if sequence != '':
        coverFiles.append((inPath+sequence+'/'+sequence+'_CoverDepthMatrix.txt'))


output = "replicon__gene"
outCover = outPath + refName + "_CoverMatrix.csv"
outDepth = outPath + refName + "_DepthMatrix.csv"

geneList = []
coverList = []
depthList = []
first_file = True
for file in coverFiles:
    (prefix, name, ext) = splitPath(file)        
    output = output + "," + name[:-17]
    coverFile = open(file)
    count = 0
    for line in coverFile:
        if line.find('\n') != -1:
            line = line[:-1]
        cover = line.split(',')
        if first_file == True:
            geneList.append(cover[0])
            coverList.append([])
            depthList.append([])
        coverList[count].append(cover[1])
        depthList[count].append(cover[2])
        count += 1
    coverFile.close()
    first_file = False

output = output + "\n"
outCoverMatrix = open(outCover, "w")
outDepthMatrix = open(outDepth, "w")
outCoverMatrix.write(output)
outDepthMatrix.write(output)

for gene in range(len(geneList)):
    outputCover = geneList[gene]
    outputDepth = geneList[gene]
    for j in range(len(coverList[gene])):
        outputCover = outputCover + "," + coverList[gene][j]
    for j in range(len(depthList[gene])):
        outputDepth = outputDepth + "," + depthList[gene][j]
    outputCover += "\n"
    outputDepth += "\n"
    outCoverMatrix.write(outputCover)
    outDepthMatrix.write(outputDepth)
    
outCoverMatrix.close()
outDepthMatrix.close()