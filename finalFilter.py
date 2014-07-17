'''
finalFilter.py

removes heterozygote calls fron a vcf, also keeping count of the het SNPs removed. 
output to _q30.vcf file and het count file

example: 
python finalFilter.py <raw>.vfc <q30>.vcf <outHetFile>

Created:	24/01/2013
Modified:	17/07/2014
author: David Edwards
'''
import sys
from pipe_utils import (splitPath)

inFile = sys.argv[1]
outFile = sys.argv[2]
(prefix,middle,ext) = splitPath(outFile)

outHetFile = sys.argv[3] + '/' + middle[:-3] + "het.txt"
vcfIn = open(inFile)
vcfOut = open(outFile, "w")
hetOut = open(outHetFile, "w")
hetCount = 0

for line in vcfIn:
    if line.startswith("#") == True:
        vcfOut.write(line)
    else:
        element = line.split("\t")        
        if (element[7].find("AF1=1") != -1 or element[-1].startswith("0") != True) and element[4].find(",") == -1:
            vcfOut.write(line)    
        elif element[7].startswith("IND") != True:
            hetCount += 1

hetOut.write(str(hetCount))
hetOut.close()
vcfOut.close()
vcfIn.close()