'''
collateRepAlleleMatrix.py

takes the SNP allele matrix for each isolate and generates the full allele matrix

outputs matrix to user-defined file

example:
python collateRepAlleleMatrix.py <input_dir> <replicon> <output>
Created:	27/2/2013
Modified:	
author: David Edwards
'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, glob
from pipe_utils import splitPath
input_dir = sys.argv[1]
replicon = sys.argv[2]
output_file = sys.argv[3]

snpList = open(snpList_name)
header = 'Pos,Ref'
(prefix1, name1, ext1) = splitPath(snpList_name)
(prefix2, name2, ext2) = splitPath(outputFile_name)

consensi = prefix1 + '/*_cns.fq'
stats = prefix2 + '/' + name2[:-8] + '_RepStats.txt'

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

# call bases for each SNP from consensus file of each strain
for file in glob.glob(consensi):
    (prefix, name, ext) = splitPath(file)
    print file 
    name = name[:-4]   
    statsFile = open(stats)
    for sample in statsFile:
        if sample.startswith("Isolate") != True:
            splitSample = sample.split("\t")
            if (name == splitSample[0]) and (splitSample[-1].startswith("f") != True):
                header = header + ',' + name
                consensus = SeqIO.parse(file, "fastq")
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
    statsFile.close()

header = header + "\n"
output = ""
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
