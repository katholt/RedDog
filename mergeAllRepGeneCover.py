'''
mergeAllRepGeneCover.py

Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)

Takes a coverage and average depth output divided by replicon__gene, and named 
after each new isolate, and converts them to an entry on the merge target 
coverage matrix and average depth matrix, divided by replicon__gene and isolate.

example:
python mergeAllGeneCover.py <inputDir> <outputDir> refName sequence_list_string

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

output = ""
CoverMatrixName = outPath + refName + "_CoverMatrix.csv"
DepthMatrixName = outPath + refName + "_DepthMatrix.csv"
geneList = []
coverList = []
depthList = []

inCoverMatrix = open(CoverMatrixName, "r")
inDepthMatrix = open(DepthMatrixName, "r")

second_line = True
for line in inCoverMatrix:
    if line.startswith('replicon__gene'):
        output = line[:-1]
    else:
        line = line[:-1]
        cover = line.split(',')
        geneList.append(cover[0])
        for counter in range(1, len(cover)):
            if second_line == True:
                coverList.append([])
            coverList[counter-1].append(cover[counter])
        second_line = False

second_line = True
for line in inDepthMatrix:
    if line.startswith('replicon__gene')!=True:
        line = line[:-1]
        depth = line.split(',')
        for counter in range(1, len(depth)):
            if second_line == True:
                depthList.append([])
            depthList[counter-1].append(depth[counter])
        second_line = False

counter = len(depthList)
for file in coverFiles:    
    (prefix, name, ext) = splitPath(file)        
    output = output + "," + name[:-17]
    coverFile = open(file)
    count = 0
    for line in coverFile:
        line = line[:-1]
        cover = line.split(',')
        if count == 0:
            coverList.append([])
            depthList.append([])
        coverList[counter].append(cover[1])
        depthList[counter].append(cover[2])
        count += 1
    coverFile.close()
    counter += 1

output = output + "\n"
outCoverMatrix = open(CoverMatrixName, "w")
outDepthMatrix = open(DepthMatrixName, "w")
outCoverMatrix.write(output)
outDepthMatrix.write(output)
for gene in range(len(geneList)):
    outputCover = geneList[gene]
    outputDepth = geneList[gene]
    for j in range(len(coverList)):
        outputCover = outputCover + "," + coverList[j][gene]
    for j in range(len(depthList)):
        outputDepth = outputDepth + "," + depthList[j][gene]
    outputCover += "\n"
    outputDepth += "\n"
    outCoverMatrix.write(outputCover)
    outDepthMatrix.write(outputDepth)    
outCoverMatrix.close()
outDepthMatrix.close()
