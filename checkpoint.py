'''
checkpoint.py

creates or appends the stage name to "checkpoint.txt" in temp folder of pipeline
indicates that the stage has successfully completed.

Note: The pipe does not use the information in "checkpoint.txt"; this is provided
for trouble-shooting pipeline problems (See manual).

example:
python checkpoint.py <temp_directory> <stage> 

Created:    14/08/2014
Modified:   
author: David Edwards
'''

import os, sys
output_prefix = sys.argv[1]
stage_string = sys.argv[2]

cp_file_name = output_prefix + "checkpoint.txt"
if not os.path.exists(cp_file_name):
	cp_file = open(cp_file_name, "w")
else:
	cp_file = open(cp_file_name, "a")
cp_file.write(stage_string+'\n')
cp_file.close()
