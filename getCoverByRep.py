'''
getCoverByRep.py

Calculates the coverage and depth by replicon for a mapped set of reads. 
Reports, in order:
replicon name, length, depth (average), cover (%)

Writes output to a user-defined file.

Examples:
python getCoverByRep.py <reference.fasta> <coverage file> <output>

python getCoverByRep.py NC_123456.fasta sampleX_cover.txt sampleX_replicon_cover.txt

Note: coverage file needs to be pre-generated pipeup by samtools (only first four columns needed)

Created:	16/01/2013
Modified:	10/09/2013
author: David Edwards
'''
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pipe_utils import get_key

reference = sys.argv[1]
coverFile = open(sys.argv[2])

chromosomes = []
output = ""

for record in SeqIO.parse(reference, "fasta"):
	chromosomes.append([record.id, len(record), 0, 0])

largest_key = 0
for chromosome in chromosomes:
	key = get_key(chromosome[0])
	if key > largest_key:
		largest_key = key

chromosome_index = []

for i in range((largest_key)):
	chromosome_index.append([])

count = 0
for chromosome in chromosomes:
	key = get_key(chromosome[0])
	chromosome_index[(key-1)].append([chromosome[0],count])
	count += 1

previous_replicon = ""
for line in coverFile:
	cover = line.split()
	if int(cover[3]) != 0:
		if cover[0] != previous_replicon:
		    key = get_key(cover[0])
		    previous_replicon = cover[0]
		for item in chromosome_index[(key-1)]:
		    if cover[0]==item[0]:
		        index=item[1]
		        chromosomes[index][2] += 1
		        chromosomes[index][3] += int(cover[3])

coverFile.close()

for chromosome in range(0, len(chromosomes)):
	lineOut = chromosomes[chromosome][0] + "\t" + str(chromosomes[chromosome][1]) + "\t"
	if chromosomes[chromosome][2] != 0:
		lineOut += str(chromosomes[chromosome][3]/float(chromosomes[chromosome][2])) + "\t"
	else:
		lineOut += "0.0\t"
	lineOut += str(chromosomes[chromosome][2]/float(chromosomes[chromosome][1]) * 100) + "\n"
	output += lineOut

outFile = open(sys.argv[3], "w")
outFile.write(output)
outFile.close()
