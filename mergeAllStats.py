#!/bin/env python
'''
mergeAllStats.py

Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)

merges the general statistics for two separate runs of the pipe.

example:
python mergeAllStats.py <new_AllStats.txt> <mergeDirectory> 

Created:	23/10/2012
Modified:	29/10/2013 for mergeAllStats (from mergeStats)
author: David Edwards
'''
import sys
from pipe_utils import splitPath

inFileName = sys.argv[1]
(inPrefix, inName, inExt) = splitPath(inFileName)
mergeFileName = sys.argv[2] + inName + inExt
outfile_allStats_name_user = sys.argv[2] + inName + "_user" + inExt
output = ""

# Combine the two sets
# First get the merge run list
mergeFile = open(mergeFileName, "r")
for line in mergeFile:
    output += line
mergeFile.close()

# then add the new run AllStats.txt leaving off the header
inFile = open(inFileName)
for line in inFile:
    items = line.split()
    if items[0] != "Isolate":
        output += line
inFile.close()

mergeFile = open(mergeFileName, "w")
mergeFile.writelines(output)
mergeFile.close()

out_lines = []
infile_allStats = open(mergeFileName,"r")
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
