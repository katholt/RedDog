'''
deriveRepAlleleMatrix.py

takes the SNP list for a replicon and the consensus sequence 
of an Isolate and generates a colummn in the allele matrix
with regard to the reference for that Isolate

outputs matrix to user-defined file

example:
python deriveRepAlleleMatrix.py <SNPList> <output> <reference.fa> <replicon> <isolate_consensus_seq> <replicon_RepStats.txt>

Created:	26/2/2014
Modified:	
author: David Edwards
'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, glob
from pipe_utils import splitPath
snpList_name = sys.argv[1]
outputFile_name = sys.argv[2]
reference_name = sys.argv[3]
replicon = sys.argv[4]
consensus_in = sys.argv[5]
stats = sys.argv[6]

snpList = open(snpList_name)
header = 'Pos,Ref'

#make a SNP list
SNP = []
for line in snpList:
    splitLine = line.split()
    SNP.append([splitLine[0]])
snpList.close()

#populate reference calls from reference fasta

references = SeqIO.parse(reference_name, "fasta")
for reference in references:
    if reference.name == replicon:
        for call in range(len(SNP)):
            base = reference.seq[int(SNP[call][0])-1]
            SNP[call].append(base.upper())

# call bases for each SNP from consensus file of the strain
(prefix, name, ext) = splitPath(consensus_in)
name = name[:-4]   
consensus = SeqIO.parse(consensus_in, "fastq")
statsFile = open(stats)
for sample in statsFile:
    if sample.startswith("Isolate") != True:
        splitSample = sample.split("\t")
        if (name == splitSample[0]) and (splitSample[-1].startswith("f") != True):
            header = header + ',' + name
            for record in consensus:
                if record.name == replicon:                        #
                    for call in range(len(SNP)):
                        if int(SNP[call][0]) <= len(record.seq):
                            base = record.seq[int(SNP[call][0])-1]
                            qualityCall = record.letter_annotations["phred_quality"][int(SNP[call][0])-1]
                            if qualityCall >= 20 and (base in ['G', 'C', 'A', 'T']):
                                SNP[call].append(base.upper())
                            else:
                                SNP[call].append('-')
                        else:
                            SNP[call].append('-')
        if (name == splitSample[0]) and (splitSample[-1].startswith("f") == True):
            header = "fail"
statsFile.close()

header = header + "\n"
output = ""
if header.startswith("fail") != True:
    for i in range(len(SNP)):
        for j in range(len(SNP[i])):
            if j < len(SNP[i]) -1:
                output = output + SNP[i][j] + ","
            else:
                output = output + SNP[i][j] + "\n"

outputFile = open(outputFile_name, "w")
outputFile.write(header)
outputFile.write(output)
outputFile.close()
