'''
collateAllStats.py
for mapping_pipe.py


collates the general statistics for set of reads that have gone through
the pipeline during a pangenome run

example:
python collateAllStats.py reference reference_name <example>_rep_cover.txt

Created:	29042013
Modified:	05052013
author: David Edwards
'''
import sys, glob
from pipe_utils import (splitPath)

reference_name = sys.argv[1]
example_rep_cover_File = sys.argv[2]

example_rep_cover = open(example_rep_cover_File)
replicon_names =[]
for line in example_rep_cover:
    entry = line.split()
    replicon_names.append(entry[0])

(prefix, middle, ext) = splitPath(example_rep_cover_File)
outfile_allStats_name = prefix[:-4] + reference_name + '_AllStats.txt'
outfile_allStats_name_user = prefix[:-4] + reference_name + '_AllStats_user.txt'

output_allStats = ''
header_allStats = 'Isolate\t'
for replicon in replicon_names:
    header_allStats += 'Cover%_' + replicon + '\t'
for replicon in replicon_names:
    header_allStats += 'Depth_' + replicon + '\t'
for replicon in replicon_names:
    header_allStats += 'Mapped%_' + replicon + '\t'
header_allStats += 'Mapped%_Total\tTotal_Reads\tInsert_Mean\tInsert_StDev\tLength_Max\tBase_Qual_Mean\tBase_Qual_StDev\tA_%\tT_%\tG_%\tC_%\tN_%\n'

inFiles = prefix + '/' + '*_AllStats.txt'
for infile in glob.glob(inFiles):
    allStatsFile = open(infile)
    output_allStats += allStatsFile.read()
    allStatsFile.close()

outfile_allStats = open(outfile_allStats_name,"w")
outfile_allStats.write(header_allStats)
outfile_allStats.writelines(output_allStats)
outfile_allStats.close()

out_lines = []
infile_allStats = open(outfile_allStats_name,"r")
test_line = infile_allStats.readline()
entries = test_line.split()
for i in range(0, len(entries)):
	out_lines.append([entries[i]])
for line in infile_allStats:
	entries = line.split()
	for i in range(0, len(entries)):
		out_lines[i].append(entries[i])
infile_allStats.close()
output = ""
for i in range(0, len(out_lines)):
	for j in range(0, len(out_lines[i])-1):
		output += out_lines[i][j] + "\t"
	output += out_lines[i][-1] + "\n"
outfile_allStats_user = open(outfile_allStats_name_user,"w")
outfile_allStats_user.write(output)
outfile_allStats_user.close()
