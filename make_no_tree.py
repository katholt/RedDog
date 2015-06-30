'''
make_no_tree.py

Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)

creates an empty tree file when no tree is required...
i.e. this is the 'non-functional' version of makeTree

example:
python make_no_tree.py <output_directory> <input_file> <out_name> 

Created:    30/06/2015
Modified:   
author: David Edwards
'''
import os, sys
out_directory = sys.argv[1]
input_file = sys.argv[2]
out_name = sys.argv[3]

if output_directory[-1] != "/":
	output_directory += "/"
output_file = out_directory + "RAxML_bestTree." + out_name
output = open(output_file, "w")

output.write("No Tree\n")
output.close()

