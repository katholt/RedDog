'''
generate_sequence_list.py

Run if there is no 'sequence_list.txt' available.
Requires the '*_AllStats_user.txt' file be available in the 
output folder.

example:
python generate_sequence_list.py <RedDog_output_folder> 

Created:	12/03/2014
Modified:
author: David Edwards
'''
import sys, glob
from pipe_utils import splitPath

output_folder_name = sys.argv[1]
if output_folder_name[-1] != '/':
    output_folder_name += '/'
user_file_name = output_folder_name + '*_AllStats_user.txt'
file_names = []
file_names = glob.glob(user_file_name)
if len(file_names) > 1:
	print "Too many AllStats_user.txt file found in " + output_folder
	sys.exit()
elif len(file_names) < 1:
	print "No AllStats_user.txt file found in " + output_folder
	sys.exit()

user_file = open(file_names[0], 'r')
sequence_list = user_file.readline()
user_file.close()
sequences = sequence_list.split()
output = ''
print sequences
for item in range(1, len(sequences)):
	print sequences[item]
	output += sequences[item] + '\n'

output_folder_name += 'sequence_list.txt'
output_file = open(output_folder_name, 'w')
output_file.write(output)
output_file.close()
