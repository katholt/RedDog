'''
dervieStats.py

collates the general statistics for set of reads that have gone through
the pipeline

example:
python deriveStats.py ref.fasta <vcf>.txt name cover_fail depth_fail mapped_fail output_stats.txt

Created:	03/03/2012
Modified:	30/06/2013
author: David Edwards
'''
import sys
from pipe_utils import splitPath
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

reference = sys.argv[1]
referenceSeqRec = SeqIO.read(reference, "fasta") 
refLength = float(len(str(referenceSeqRec.seq)))

name = sys.argv[3]
output = name + "\t"
(prefix, middle, ext) = splitPath(sys.argv[2])
name = prefix + "/" + name

coverFile = open(name + "_ave_cover.txt")
coverage = coverFile.readline()
coverList = coverage.split()
aveReadDepth = coverList[1]
totalBases = float(coverList[2])
minDepthBases = float(coverList[3])

coverRef = 100 * totalBases/refLength
if totalBases == 0:
	minCover = 0
else:
	minCover = 100 * minDepthBases/totalBases
output = output + str(coverRef) +"\t"+ aveReadDepth +"\t"+ str(minCover) +"\t"

mapped = open(name + "_bam.txt")
mappedList = mapped.readlines()
totalReadsList = mappedList[5].split()
totalReads = float(totalReadsList[2])
mappedReadsList = mappedList[6].split()
mappedReads = float(mappedReadsList[2])
mappedRef = 100 * mappedReads/totalReads
output = output + str(int(totalReads)) +"\t"+ str(mappedRef) +"\t"

vcfFile = open(sys.argv[2])
hetFile = open(name + "_het.txt")
vcf = vcfFile.readline()
het = hetFile.read()
vcfList = vcf.split()
output = output + str(int(vcfList[0])) +"\t"+ str(int(het)) +"\t"+ str(int(vcfList[2]))

# decide if a sample (strain) is a fail or pass based on these three stats
coverFail = int(sys.argv[4])
depthFail = int(sys.argv[5])
mappedFail = int(sys.argv[6])
if (coverRef < coverFail):
	output = output + "\t" + "f"
elif (float(aveReadDepth) < depthFail):
	output = output + "\t" + "f"
elif (mappedRef < mappedFail):
	output = output + "\t" + "f"	
output = output +"\n"
outFile = open(sys.argv[7], "w")
outFile.write(output)
outFile.close()


