'''
collateRepAlleleMatrix.py

Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)

takes the SNP allele matrix entry for each isolate and generates the full allele matrix

outputs matrix to user-defined file

example:
python collateRepAlleleMatrix.py <temp_dir> <output> sequences_string rep_name 

Created:	27/2/2014
Modified:	21/7/2014
author: David Edwards
'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, glob
from rubra.utils import splitPath

prefix = sys.argv[1]
output_file = sys.argv[2]
sequences_string = sys.argv[3]
rep_name = sys.argv[4]

input_files = []
sequences = sequences_string.split(',')
for sequence in sequences:
    if sequence != '':
        input_files.append((prefix+sequence+'/deriveRepAlleleMartix/'+rep_name+'_'+sequence+'_alleles.txt'))

header = ''
SNPmatrix = []

for entry_file in input_files:
    matrix_entry_file = open(entry_file, 'r')
    matrix_entry = matrix_entry_file.readline()
    if matrix_entry.startswith('fail') != True:
        if header == '':
            header = matrix_entry[:-1]
            matrix_entrys = matrix_entry_file.readlines()
            for line in matrix_entrys:
                SNPmatrix.append([line[:-1]])

        else:
            new_isolate = matrix_entry.split(',')
            header += ',' + new_isolate[2][:-1]
            matrix_entrys = matrix_entry_file.readlines()
            count = 0
            for line in matrix_entrys:
                entry = line.split(',')
                new_string = SNPmatrix[count][0]
                new_string += "," + entry[2][:-1]
                SNPmatrix[count][0] = new_string
                count += 1
    matrix_entry_file.close()

if header == '':
    header = 'Pos,Ref'
header = header + "\n"
output = ""
for i in range(len(SNPmatrix)):
    output += SNPmatrix[i][0] + "\n"

outputFile = open(output_file, "w")
outputFile.write(header)
outputFile.write(output)
outputFile.close()
