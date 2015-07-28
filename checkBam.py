#!/bin/env python
'''
checkBam.py

Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)

simple check of the constructed BAM:
For PE reads, if the BAM is not larger or equal to the size of
one of the read sets (_1), the BAM is deleted and an error code returned
For SE and IT reads, if the BAM is less than half the size of the read set,
the BAM is deleted and an error code returned - this should stop the pipeline!

otherwise, nothing happens...

example:
python checkBam.py <bam> <read_type> <sequence_mapped> 

[only need one of the pair for PE reads]

Created:    24/07/2015
Modified:   
author: David Edwards
'''
import os, sys
bam_to_check = sys.argv[1]
read_type = sys.argv[2]
sequence_mapped = sys.argv[3]
try:
	bam_info = os.stat(bam_to_check)
except:
	print "BAM file does not exist"
	sys.exit(1)
seq_info = os.stat(sequence_mapped)
if read_type == 'PE':
	if bam_info.st_size < seq_info.st_size:
		os.remove(bam_to_check)
		print "BAM too small: deleted"
		sys.exit(1)
else:
	if bam_info.st_size < (seq_info.st_size*0.5):
		os.remove(bam_to_check)
		print "BAM too small: deleted"
		sys.exit(1)
