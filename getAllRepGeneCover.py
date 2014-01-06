'''
getAllRepGeneCover.py

Takes the genbank file for the reference sequence and all the coverage files 
(from bamtools), and converts then to a coverage matrix and average depth matrix,
divided by replicon__gene and strain.

example:
python getAllGeneCover.py <inputDir> <outputDir> ref.genbank

Created:	17/5/2013
Modified:   26/5/2013
author: David Edwards
'''
import sys, glob
from pipe_utils import splitPath, get_key
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC

inPath =sys.argv[1]
outPath = sys.argv[2]
genbankName = sys.argv[3]

(refPrefix, refName, refExt) = splitPath(genbankName)
output = "replicon__gene"
coverFileName = inPath + "temp/*_coverage.txt"
outCover = outPath + refName + "_CoverMatrix.csv"
outDepth = outPath + refName + "_DepthMatrix.csv"

geneList = []
geneCoverList = []
geneDepthList = []
totalBases = 0
repliconList =[]

records = SeqIO.parse(genbankName, "genbank")
for record in records:
    feature_count = 0
    for f in record.features:
        if f.type == "CDS":
            feature_count += 1
            start = f.location.nofuzzy_start
            stop = f.location.nofuzzy_end
            geneList.append([record.name, start, stop, 0, 0])
            sysid = record.name+";"+str(start)+'-'+str(stop)
            f.qualifiers['sysid'] = [sysid]
            if 'locus_tag' in f.qualifiers:
                locus_tag = f.qualifiers['locus_tag'][0]
            else:
                #if the locus_tag is missing from the genbank record make up a tag
                locus_tag = "tag_" + str(start)+'-'+str(stop)
            geneCoverList.append([record.name+'__'+locus_tag])
            geneDepthList.append([record.name+'__'+locus_tag])
    repliconList.append([record.name,len(record),feature_count])

largest_key = 0
for replicon in repliconList:
    key = get_key(replicon[0])
    if key > largest_key:
        largest_key = key
replicon_index = []
for i in range((largest_key)):
    replicon_index.append([])
count = 0
for replicon in repliconList:
    key = get_key(replicon[0])
    replicon_index[(key-1)].append([replicon[0],count])
    count += 1

geneSlice = []
sliceSize = []
#slice the geneList by coordinates and replicon
count = 0
for replicon in repliconList:
    geneSlice.append([])
    if replicon[2] > 0:
        sliceSize.append(replicon[1] / replicon[2] + 1)
    else:
        sliceSize.append(2)
    for slice in range(replicon[2]+1):
        geneSlice[count].append([])
    count+=1

geneCount = 0
for gene in geneList:
    key = get_key(gene[0])
    for item in replicon_index[(key-1)]:
            if gene[0]==item[0]:
                index=item[1]
    slice1 = geneList[geneCount][1]/sliceSize[index]
    slice2 = geneList[geneCount][2]/sliceSize[index]
    geneSlice[index][slice1].append([geneList[geneCount][0], geneList[geneCount][1], geneList[geneCount][2], geneCount])
    while slice1 < slice2:
        slice1 += 1
        geneSlice[index][slice1].append([geneList[geneCount][0], geneList[geneCount][1], geneList[geneCount][2], geneCount])
    geneCount += 1


for file in glob.glob(coverFileName):
    (prefix, name, ext) = splitPath(file)        
    output = output + "," + name[:-9]
    #calculate the cover and depth for slice of genes
    coverFile = open(file)
    for line in coverFile:
        cover = line.split()
        key = get_key(cover[0])
        for item in replicon_index[(key-1)]:
            if cover[0]==item[0]:
                index=item[1]
        if int(cover[2]) > 0:
            slice = (int(cover[1]) + 1)/sliceSize[index]
            geneCount = 0
            if geneSlice[index][slice] != []:
                for gene in geneSlice[index][slice]:
                    if geneSlice[index][slice][geneCount][0]==cover[0]:
                        if int(cover[1]) + 1 >= geneSlice[index][slice][geneCount][1]:
                            if int(cover[1]) + 1 <= geneSlice[index][slice][geneCount][2]:
                                geneList[geneSlice[index][slice][geneCount][3]][3] += 1
                                geneList[geneSlice[index][slice][geneCount][3]][4] += int(cover[2])                        
                    geneCount += 1

    coverFile.close()
    geneCount = 0
    for gene in geneList:
        coverage = float(geneList[geneCount][3]) * 100.0 / float(geneList[geneCount][2] - geneList[geneCount][1] + 1)
        if geneList[geneCount][3]!=0:
            depth = float(geneList[geneCount][4]) / float(geneList[geneCount][3])
        else:
            depth = 0.0
        geneCoverList[geneCount].append(coverage)
        geneDepthList[geneCount].append(depth)
        geneList[geneCount][3] = 0
        geneList[geneCount][4] = 0                        
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
