#!/bin/env python
'''
makeDir.py

Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)

creates all the output and temp folders required by the pipeline

example:
python makeDir.py <output_directory> <sequences_string> 

Created:    14/07/2014
Modified:   17/07/2014
author: David Edwards
'''
import os, sys
output_prefix = sys.argv[1]

sequences_string = sys.argv[2]
sequences = sequences_string.split(',')

seq_stages = ["callRepSNPs", "q30VarFilter/hets", "getVCFStats", "deriveRepStats", "deriveRepAlleleMartix"]
temp_stages = ["getRepSNPList", "collateRepAlleleMatrix"]

dir_names = []

dir_names.append(output_prefix + "bam/")
dir_names.append(output_prefix + "vcf/")
dir_names.append(output_prefix + "temp/success/")

for temp_stage in temp_stages:
	dir_names.append(output_prefix + "temp/" + temp_stage + "/")

for sequence in sequences:
	if sequence != "":
		for seq_stage in seq_stages:
			dir_names.append(output_prefix + "temp/" + sequence + "/" + seq_stage + "/")

for dir_name in dir_names:
	if not os.path.isfile(dir_name) and not os.path.isdir(dir_name):
		os.makedirs(dir_name)

