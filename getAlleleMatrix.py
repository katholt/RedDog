'''
getAlleleMatrix.py

takes the SNP list and the consensus sequences of all ingroups 
and outgroups and generates an allele matrix with regard to 
the reference

outputs matrix to user-defined file

example:
python getSNPList.py <SNPList> <output> <reference.fa>

Created:	15/5/2011
Modified:	22/2/2013
author: David Edwards
'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, glob
from pipe_utils import splitPath

snpList = open(sys.argv[1])
header = 'Pos,Ref'
(prefix1, name1, ext1) = splitPath(sys.argv[1])
(prefix2, name2, ext2) = splitPath(sys.argv[2])

consensi = prefix1 + '/*_cns.fq'
print consensi
stats = prefix2 + '/' + name2[:-8] + '_stats.tab'
print stats
#make a SNP list
SNP = []
for line in snpList:
    splitLine = line.split()
    SNP.append([splitLine[0]])
snpList.close()
#populate reference calls from reference fasta
reference = SeqIO.parse(sys.argv[3], "fasta").next()
for call in range(len(SNP)):
    base = reference.seq[int(SNP[call][0])-1]
    SNP[call].append(base.upper())
# call bases for each SNP from consensus file of each strain
for file in glob.glob(consensi):
    (prefix, name, ext) = splitPath(file)
    statsFile = open(stats)
    for sample in statsFile:
        if sample.startswith("Name") != True:
            splitSample = sample.split("\t")
            if (name[:-4] == splitSample[0]) and (splitSample[-1].startswith("f") != True):
                header = header + ',' + name[:-4]
                consensus = SeqIO.parse(file, "fastq").next()
                for call in range(len(SNP)):
                    test = "GACTgact"
                    if int(SNP[call][0]) <= len(consensus.seq):
                        base = consensus.seq[int(SNP[call][0])-1]
                        qualityCall = consensus.letter_annotations["phred_quality"][int(SNP[call][0])-1]
                        if qualityCall >= 20 and test.find(base) != -1:
                            SNP[call].append(base.upper())
                        else:
                            SNP[call].append('-')
                    else:
                        SNP[call].append('-')
    statsFile.close()
header = header + "\n"
outputFile = open(sys.argv[2], "w")
outputFile.write(header)
for i in range(len(SNP)):
    output = ""
    for j in range(len(SNP[i])):
        if j < len(SNP[i]) -1:
            output = output + SNP[i][j] + ","
        else:
            output = output + SNP[i][j] + "\n"
    outputFile.write(output)
outputFile.close()