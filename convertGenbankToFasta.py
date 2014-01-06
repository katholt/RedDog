'''
convertGenbankToFasta.py

Converts a genbank file to a fasta file 
  - produces a sequence for each replicon
  - uses the sequence name to generate the id for the seqeunce
  - writes output to a user-defined fasta file.

Examples:
python convertGenbankToFasta.py <genbank> <output>

Created:	07052013
Modified:	
author: David Edwards
'''
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

input_handle = open(sys.argv[1], "rU")
output_handle = open(sys.argv[2], "w")
sequences_out = []
sequences = SeqIO.parse(input_handle, "genbank")
for sequence in sequences:
	new_seq = SeqRecord(sequence.seq)
	new_seq.id = sequence.name
	new_seq.description = sequence.description
	sequences_out.append(new_seq)
SeqIO.write(sequences_out, output_handle, "fasta")
 
output_handle.close()
input_handle.close()
