#!/bin/env python
'''
make_no_tree.py

Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)

creates an empty tree file when no tree is required...
i.e. this is the 'non-functional' version of makeTree

example:
python make_no_tree.py <output_directory> <input_file> <out_name> 

Created:    30/06/2015
Modified:   27/10/2015 (change back to FastTree)
author: David Edwards
'''
import os, sys
input_file = sys.argv[1]
output_file = sys.argv[2]

output = open(output_file, "w")
output.write("No Tree\n")
output.close()

