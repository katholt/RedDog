#!/usr/bin/env python
'''
deriveAllRepGeneCover.py

Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)

Takes the genbank file for the reference sequence and the coverage file for an isolate 
(from bamtools), and converts them to a coverage and average depth output,
divided by replicon__gene, and named after the isolate for further collation.

example:
python deriveAllGeneCover.py <outputDir> ref.genbank <coverfile>

Created:  27/05/2013
Last Modified: 17/02/2015 - fixed locus tag designation
author: David Edwards
'''
import sys, glob
from pipe_utils import splitPath, get_key
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC

outPath = sys.argv[1]
genbankName = sys.argv[2]
coverFileName = sys.argv[3]

(prefix, name, ext) = splitPath(coverFileName)        
outCover = outPath + name[:-9] + "/" + name[:-9] + "_CoverDepthMatrix.txt"

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
#        if f.type == "misc_feature":
            feature_count += 1
            start = f.location.nofuzzy_start
            stop = f.location.nofuzzy_end
            geneList.append([record.name, start, stop, 0, 0])
            sysid = record.name+";"+str(start+1)+'-'+str(stop+1)
            f.qualifiers['sysid'] = [sysid]
            if 'locus_tag' in f.qualifiers:
                locus_tag = f.qualifiers['locus_tag'][0]
            else:
                #if the locus_tag is missing from the genbank record make up a tag
                locus_tag = "tag_" + str(start+1)+'-'+str(stop+1)
            geneCoverList.append([record.name+'__'+locus_tag])
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
        sliceSize.append(replicon[1] + 1)
    if replicon[2] > 0:
        for slice in range(replicon[2]+1):
            geneSlice[count].append([])
    else:
        geneSlice[count].append([])
        geneSlice[count].append([])
    count+=1

geneCount = 0
previous_replicon = ""
for gene in geneList:
    if gene[0] != previous_replicon:
        key = get_key(gene[0])
        previous_replicon = gene[0]
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

#calculate the cover and depth for slice of genes
previous_replicon = ""
coverFile = open(coverFileName)
for line in coverFile:
    cover = line.split()
    if int(cover[3]) > 0:
        if cover[0] != previous_replicon:
            key = get_key(cover[0])
            previous_replicon = cover[0]
            for item in replicon_index[(key-1)]:
                    if cover[0]==item[0]:
                        index=item[1]
        slice = (int(cover[1]) + 1)/sliceSize[index]
        geneCount = 0
        if geneSlice[index][slice] != []:
            for gene in geneSlice[index][slice]:
                if geneSlice[index][slice][geneCount][0]==cover[0]:
                    if int(cover[1]) + 1 >= geneSlice[index][slice][geneCount][1]:
                        if int(cover[1]) + 1 <= geneSlice[index][slice][geneCount][2]:
                            geneList[geneSlice[index][slice][geneCount][3]][3] += 1
                            geneList[geneSlice[index][slice][geneCount][3]][4] += int(cover[3])                        
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
    geneCoverList[geneCount].append(depth)
    geneCount += 1

outCoverMatrix = open(outCover, "w")
for gene in range(len(geneList)):
    outputCover = geneCoverList[gene][0]+ "," + str(geneCoverList[gene][1]) + "," + str(geneCoverList[gene][2]) + "\n"
    outCoverMatrix.write(outputCover)
outCoverMatrix.close()

