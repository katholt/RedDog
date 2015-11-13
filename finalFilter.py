#!/bin/env python
'''
finalFilter.py

Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)

Filters out snp calls with qual less than 30 and removes heterozygote calls fron a vcf, 
also keeping count of the het SNPs removed. 
output to _q30.vcf file and het count file

If HetsVCF is set to 'True', the heterozygous SNP calls will be written to a VCF

example: 
python finalFilter.py <raw>.vfc <q30>.vcf <outHetFile> <HetsVCF>

Created:	24/01/2013
Modified:	13/11/2015
author: David Edwards
'''
import sys
from pipe_utils import (splitPath)

inFile = sys.argv[1]
outFile = sys.argv[2]
(prefix,middle,ext) = splitPath(outFile)

outHetFile = sys.argv[3] + '/' + middle[:-3] + "het.txt"
if sys.argv[4] == "True":
    HetsVCF = True
else:
    HetsVCF = False
vcfIn = open(inFile)
vcfOut = open(outFile, "w")
hetOut = open(outHetFile, "w")
if HetsVCF:
    outHetPosFile = prefix + '/' + middle[:-3] + "het.vcf"
    hetVcfOut = open(outHetPosFile,"w")
hetCount = 0

for line in vcfIn:
    if line.startswith("#") == True:
        vcfOut.write(line)
        if HetsVCF:
        	hetVcfOut.write(line)
    else:
        element = line.split("\t")
        if float(element[5]) >= 30:       
            if (element[7].find("AF1=1") != -1 or element[-1].startswith("0") != True) and element[4].find(",") == -1:
                vcfOut.write(line)    
            elif element[7].startswith("IND") != True:
                hetCount += 1
                if HetsVCF:
                    hetVcfOut.write(line)

hetOut.write(str(hetCount))
hetOut.close()
vcfOut.close()
vcfIn.close()
if HetsVCF:
	hetVcfOut.close()
