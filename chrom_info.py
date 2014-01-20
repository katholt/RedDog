import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def chromInfoFasta(refFile):
    # extract the chromosome name and length from the fasta reference
    chroms = []
    for record in SeqIO.parse(refFile, "fasta"):
        name = record.name
        if name.find('.') != -1:
            temp_name = name.split('.')
            name = temp_name[0]
        chroms.append((name, len(record)))
    return chroms

def chromInfoGenbank(refFile):
    # extract the chromosome name and length from the genbank reference
    chroms = []
    for record in SeqIO.parse(refFile, "genbank"):
        name = record.name
        if name.find('.') != -1:
            temp_name = name.split('.')
            name = temp_name[0]
        chroms.append((name, len(record)))
    return chroms
