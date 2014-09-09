#!/usr/bin/env python
#
# Converts Mummer repeat-match output to coordinate table for parseSNPTable
#
# Author(s) - David Edwards
#
# Example command:
'''
module load python-gcc/2.7.5
python convertRepeatsToCoords.py -i repeats.txt -o coords.txt
'''
#
# Created: 20140909
# Changes:
#	 <date>   - <change>
#

import os, sys, glob
import subprocess
import string
import re
#import random
import operator
from optparse import OptionParser
#import resource

def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i", "--input", action="store", dest="input", help="Input: repeats file from Mummer (default none)", default="")

	parser.add_option("-o", "--output", action="store", dest="output", help="Output: coords file for parseSNPTable (default none)", default="")

	
	return parser.parse_args()


if __name__ == "__main__":

	(options, args) = main()
	
	def getRepeats(repeats, input):

		in_file = open(input)
		for line in in_file:
			if line.startswith("Long") != True and line.startswith("   Start1") != True:
				while line.startswith(" ") == True:
					line = line.lstrip(" ")
				repeat = line.split()
				size = int(repeat[2])
				start1 = int(repeat[0])
				stop1 = start1 + size
				if repeat[1].find("r") == -1:
					start2 = int(repeat[1])
					stop2 = start2 + size
				else:
					start2 = int(repeat[1].rstrip("r"))
					stop2 = start2 - size
				repeats.append([str(start1),str(stop1)])
				repeats.append([str(start2),str(stop2)])
		in_file.close()
		return(repeats)

	def repeatsOut(repeats, output):

		repeats_out = ""
		count = 0
		for repeat in repeats:
			repeats_out += (repeat[0] + "," + repeat[1] +"\n")
		output_file = open(output, "w")
		output_file.write(repeats_out)
		output_file.close()
		return
		
	### MAIN PROCESS
			
	# set up variables
	if options.input == '':
		print '\nNo input repeats file provided (-i)'
		sys.exit()
	else:
		input = options.input
	if options.output == '':
		print '\nNo output coords file provided (-o)'
		sys.exit()
	else:
		output = options.output

	repeats = []
	repeats = getRepeats(repeats, input)
	repeatsOut(repeats, output)
	
#	print ('Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
