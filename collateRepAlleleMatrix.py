'''
collateRepAlleleMatrix.py

takes the SNP allele matrix entry for each isolate and generates the full allele matrix

outputs matrix to user-defined file

example:
python collateRepAlleleMatrix.py <input> <output> length_to_remove 

Created:	27/2/2014
Modified:	18/3/2014
author: David Edwards
'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, glob
from rubra.utils import splitPath

input_file = sys.argv[1]
length_to_remove = int(sys.argv[3])*-1
(prefix, name, ext) = splitPath(input_file)
input_files = prefix + '/' + name[:length_to_remove] + '*_alleles.txt'

output_file = sys.argv[2]
header = ''
SNPmatrix = []

for entry_file in glob.glob(input_files):
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
