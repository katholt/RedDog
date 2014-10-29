'''
getVcfStats.py
based on vcfStats.py

counts number of SNPs and INDELS in _q30.vcf files
output to user-defined file

example:
python getVcfStats.py vcfFile.vcf output.txt

Created:	1/3/2012 (19/10/2011)
Modified:	18/02/2013
author: David Edwards
'''
import sys

inFile30 = sys.argv[1]
vcfFile30 = open(inFile30)
outputFile = open(sys.argv[2], "w")
snpCount30= indelCount30 = twoAltCount30 = 0
for line in vcfFile30:
    if line.startswith("#") == False:
        splitLine = line.split("\t")        
        if splitLine[7].startswith("I") == True:
            indelCount30 += 1
        else:
            snpCount30 += 1
            if splitLine[4].find(",") != -1:
                twoAltCount30 += 1

outputFile.write(str(snpCount30) + " " + str(twoAltCount30) + " " + str(indelCount30))

outputFile.close()
vcfFile30.close()
