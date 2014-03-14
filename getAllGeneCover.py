'''
getAllGeneCover.py
for pipe_AllGeneCover.py
(further modified for pipe_vda.py)
(further modified for mapping_pipe.py)

Takes the genbank file for the reference sequence and all the coverage files 
(from bamtools), and converts then to a coverage matrix and average depth matrix,
divided by gene and strain. Failed set are now excluded.

example:
python getAllGeneCover.py <inputDir> <outputDir> ref.genbank

Created:	04/5/2012
Modified:	21/3/2013
author: David Edwards
'''
import sys, glob
from pipe_utils import splitPath
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
#note: doTest is to allow for user getting output for only those samples that have passed

doTest=False
inPath =sys.argv[1]
outPath = sys.argv[2]
genbankName = sys.argv[3]

(refPrefix, refName, refExt) = splitPath(genbankName)
output = "gene"
coverFileName = inPath + "temp/*_coverage.txt"
statsFileName = outPath + refName + "_stats.tab"
outCover = outPath + refName + "_CoverMatrix.csv"
outDepth = outPath + refName + "_DepthMatrix.csv"

geneList = []
geneCoverList = []
geneDepthList = []
totalBases = 0

record = SeqIO.parse(genbankName, "genbank").next()
totalBases = len(record)
for f in record.features:
    if f.type == "CDS":
        start = f.location.nofuzzy_start
        stop = f.location.nofuzzy_end
        geneList.append([start, stop, 0, 0])
        sysid = record.name+";"+str(start)+'-'+str(stop)
        f.qualifiers['sysid'] = [sysid]
        if 'locus_tag' in f.qualifiers:
            locus_tag = f.qualifiers['locus_tag'][0]
        else:
            #if the locus_tag is missing from the genbank record make up a tag
            locus_tag = "tag_" + str(start)+'-'+str(stop)
        geneCoverList.append([locus_tag])
        geneDepthList.append([locus_tag])

#slice the geneList by coordinates
sliceSize = totalBases / 10000 + 1
geneSlice = []
for slice in range(10001):
    geneSlice.append([])

geneCount = 0
for gene in geneList:
    slice1 = geneList[geneCount][0]/sliceSize
    slice2 = geneList[geneCount][1]/sliceSize
    geneSlice[slice1].append([geneList[geneCount][0], geneList[geneCount][1], geneCount])
    while slice1 < slice2:
        slice1 += 1
        geneSlice[slice1].append([geneList[geneCount][0], geneList[geneCount][1], geneCount])
    geneCount += 1

for file in glob.glob(coverFileName):
    (prefix, name, ext) = splitPath(file)        
    statsFile = open(statsFileName)
    for sample in statsFile:
        if sample.startswith("Name") != True:
            splitSample = sample.split("\t")
            if (name[:-9] == splitSample[0]) and ((splitSample[-1].startswith("f") != True) or (doTest!=True)):
                output = output + "," + name[:-9]
                #calculate the cover and depth for slice of genes
                coverFile = open(file)
                for line in coverFile:
                    cover = line.split()    
                    if int(cover[2]) > 0:
                        slice = (int(cover[1]) + 1)/sliceSize
                        geneCount = 0
                        if geneSlice[slice] != []:
                            for gene in geneSlice[slice]:
                                if int(cover[1]) + 1 >= geneSlice[slice][geneCount][0]:
                                    if int(cover[1]) + 1 <= geneSlice[slice][geneCount][1]:
                                        geneList[geneSlice[slice][geneCount][2]][2] += 1
                                        geneList[geneSlice[slice][geneCount][2]][3] += int(cover[2])                        
                                geneCount += 1
                coverFile.close()
                geneCount = 0
                for gene in geneList:
                    cover = float(geneList[geneCount][2]) * 100.0 / float(geneList[geneCount][1] - geneList[geneCount][0] + 1)
                    if geneList[geneCount][2]!=0:
                        depth = float(geneList[geneCount][3]) / float(geneList[geneCount][2])
                    else:
                        depth = 0.0
                    geneCoverList[geneCount].append(cover)
                    geneDepthList[geneCount].append(depth)
                    geneList[geneCount][2] = 0
                    geneList[geneCount][3] = 0                        
                    geneCount += 1

output = output + "\n"
outCoverMatrix = open(outCover, "w")
outDepthMatrix = open(outDepth, "w")
outCoverMatrix.write(output)
outDepthMatrix.write(output)

for gene in range(len(geneList)):
    outputCover = geneCoverList[gene][0]
    outputDepth = geneDepthList[gene][0]
    for j in range((len(geneCoverList[gene])-1)):
        outputCover = outputCover + "," + str(geneCoverList[gene][j+1])
    for j in range((len(geneDepthList[gene])-1)):
        outputDepth = outputDepth + "," + str(geneDepthList[gene][j+1])
    outputCover += "\n"
    outputDepth += "\n"
    outCoverMatrix.write(outputCover)
    outDepthMatrix.write(outputDepth)
    
outCoverMatrix.close()
outDepthMatrix.close()
