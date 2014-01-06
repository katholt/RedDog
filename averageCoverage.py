'''
averageCoverage.py
for pipe_VariantDiscovery.py

calculates the average coverage (depth) for a mapped set of reads: reports (in order)
twice the average (as integer), the average (as a float), the number of bases with at least one mapped read and number of bases with at least a certain minimum depth of reads.  

example:
python averageCoverage.py sample_coverage.txt minDepth sample_ave_cover.txt

Created:	03/3/2012
Modified:	30/6/2013
author: David Edwards
'''
import sys
from pipe_utils import splitPath
total = 0
count = 0
minDepthCount = 0
sample = open(sys.argv[1])
minDepth = int(sys.argv[2])
for line in sample:
    data = line.split()
    if int(data[2]) != 0:
        count += 1
        total += int(data[2])
        if int(data[2]) >= minDepth:
            minDepthCount += 1
if count == 0:
	output = "0 0 0 " + str(minDepthCount)
else:
	output = str(int((total/count)*2+0.5)) +" "+ str(total/(count*1.0)) +" "+ str(count) +" "+ str(minDepthCount)
outFile = open(sys.argv[3], "w")
outFile.write(output)
outFile.close()
