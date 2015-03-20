'''
Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)
'''
#!/usr/bin/env python
#
# Converts RedDog coverage putput to heatmap data for plotTree.R
#
# Need to specify either by replicon (-i AllStats.txt) or by gene coverage (-I CoverMatrix.csv), or both
# If both are specified, coverage in the heat map matrix will reported for each replicon then its genes
#
# Replicons (-r) is only required if there is more than one replicon reported in the data set
# Reporting of replicons (and/or their genes) will be ordered as entered in -r (Comma separated, no spaces)
#
# For replicons, any isolate with read depth less than -d (default 10) will be set to 0% coverage (set to 0 to negate)
#
# For genes, the user can set a threshold coverage (-c); if the gene coverage is lower than the threshold,
# the isolate will be set to 0% coverage for that gene
#
# Author(s) - David Edwards
#
# Example command:
'''
module load python-gcc/2.7.5
python get_cover.py -l include_list.txt -i my_AllStats.txt -d 10 -I my_GeneCover.txt -c 95 -o my_heatmap.csv -r replicon3,replicon2,replicon4
'''
#
# Created: 20141021
# Changes:
#	 <date>   - <change>
#

import os, sys, glob
import subprocess
import string
import re
import operator
from optparse import OptionParser
#import resource

def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-r", "--replicons", action="store", dest="replicons", help="Replicon(s) to report, if all replicons not required - also determines reporting order (default: none)", default="")
	parser.add_option("-l", "--include_list", action="store", dest="include_list", help="List of isolates to report, otherwise all isolates included (default: none)", default="")

	parser.add_option("-i", "--input_replicon", action="store", dest="input_replicon", help="Input: coverage statistics file for replicons (required*: RedDog AllStats.txt)", default="")
	parser.add_option("-d", "--depth", action="store", dest="depth", help="Minimum required depth for replicon(s) to report coverage, otherwise reports 0% coverage for that isolate/replicon pair (default: 10)", default="10")

	parser.add_option("-I", "--input_gene", action="store", dest="input_gene", help="Input: gene cover file (required*: RedDog CoverMatrix.csv)", default="")
	parser.add_option("-c", "--coverage", action="store", dest="coverage", help="Minimum required coverage for genes to report coverage, otherwise reports 0% coverage for that isolate/gene pair (default: none, reports value in gene cover file)", default="")

	parser.add_option("-o", "--output", action="store", dest="output", help="Output: csv file for generating heatmap in plotTree.R (required)", default="")
	
	return parser.parse_args()


if __name__ == "__main__":

	(options, args) = main()
	
	def getHeatMapData(replicons_out, genes_out, replicons, input_replicon, input_gene, include_list, test_depth, test_coverage):

		include = []
		genes_out_index = []
		genes_list = []
		if include_list != "":
			try:
				includes = open(include_list)
				for line in includes:
					include.append(line.rstrip('\n'))
				includes.close()
			except:
				print "\nCould not get isolate list: check option -l"
				os._exit()

		if input_replicon != "":
			try:
				replicon_file = open(input_replicon)
				replicon_index = []
				replicon_list = []
				for line in replicon_file:
					if line.startswith('Isolate'):
						data = line.rstrip('\n').split('\t')
						count = 0
						for data_entry in data:
							if data_entry.startswith('Cover%'):
								count += 1
								replicon_name = data_entry[7:]
								if replicons != []:
									if replicon_name in replicons:
										replicon_list.append(replicon_name)
										replicon_index.append(count)
								else:
									replicon_list.append(replicon_name)
									replicon_index.append(count)
					else:
						if replicon_index == []:
							print "\nCould not find any replicons: check option -r"
							replicon_file.close()
							os._exit()						
						# order replicon_index on replicons (if not [])
						ordered_replicon_index = []
						ordered_replicon_list = []
						if replicons != [] and len(replicon_index) > 1:
							for replicon in replicons:
								if replicon in replicon_list:
									ordered_replicon_index.append(replicon_index[replicon_list.index(replicon)])
									ordered_replicon_list.append(replicon)
							replicon_index = ordered_replicon_index
							replicon_list = ordered_replicon_list
						replicons_out_string = ""
						data = line.split('\t')
						if include == [] or data[0] in include:
							replicons_out_string += data[0]
							for i in range(len(replicon_index)):
								if float(data[replicon_index[i]+count]) >= test_depth:
									replicons_out_string += (',' + replicon_list[i] + ',' + data[replicon_index[i]])
								else:
									replicons_out_string += (',' + replicon_list[i] + ',0')
						if replicons_out_string != "":
							replicons_out.append(replicons_out_string)
				replicon_file.close()
			except:
				print "\nCould not open AllStats file: check option -i " + input_replicon
				sys.exit()

		if input_gene != "":
			try:
				gene_file = open(input_gene)
				isolate_index = []
				isolate_list = []
				replicon_list = []
				for line in gene_file:
					if line.startswith('replicon__gene'):
						line = line.rstrip()
						data = line.split(',')
						for i in range(1,len(data)):
							if include == [] or data[i] in include:
								isolate_list.append(data[i])
								isolate_index.append(i)
					else:
						if isolate_index == []:
							print "\nCould not find any isolates: check isolate list file (-l)"
							gene_file.close()
							os._exit()						
						line = line.rstrip()
						data = line.split(',')
						replicon_name = data[0].split('__')
						if replicons == [] or replicon_name[0] in replicons:
							if replicon_name[0] not in replicon_list:
								replicon_list.append(replicon_name[0])
								genes_list.append([replicon_name[0], replicon_name[1]])
								for i in range(len(isolate_list)):
									test_name = replicon_name[0]+'__'+isolate_list[i]
									if float(data[isolate_index[i]]) >= test_coverage:
										genes_out.append([replicon_name[0], isolate_list[i], data[isolate_index[i]]])
										genes_out_index.append(test_name)
									else:
										genes_out.append([replicon_name[0], isolate[i], "0"])
										genes_out_index.append(test_name)
							else:
								genes_list[replicon_list.index(replicon_name[0])].append(replicon_name[1])
								for i in range(len(isolate_list)):
									test_name = replicon_name[0]+'__'+isolate_list[i]
									if float(data[isolate_index[i]]) >= test_coverage:
										genes_out[genes_out_index.index(test_name)].append(data[isolate_index[i]])
									else:
										genes_out[genes_out_index.index(test_name)].append("0")
				gene_file.close()
			except:
				print "Could not open Gene Summary file: check option -I " + input_gene
				sys.exit()
		return(replicons_out, genes_out, genes_out_index, genes_list)

	def heatMapDataOut(replicons_out, genes_out, genes_out_index, genes_list, replicons, output):
		
		output_string = ""
		header = "Isolate"
		replicon_list = []

		if replicons_out != [] and genes_out == []:
			for report in replicons_out:
				data = report.split(',')
				output_string += data[0]
				if data[1] not in replicon_list:
					header += ',' + data[1]
					replicon_list.append(data[1])
				for i in range(2, len(data), 2):
					output_string += ("," + data[i])
				output_string += "\n"

		if replicons_out != [] and genes_out != []:
			for report in replicons_out:
				data = report.split(',')
				output_string += data[0]
				for i in range(1, len(data), 2):
					if data[i] not in replicon_list:
						header += ',' + data[i]
						replicon_list.append(data[i])
						for report in genes_list:
							if report[0] == data[i]:
								for j in range(1,len(report)):
									header += (',' + report[j])
				for i in range(2, len(data), 2):
					output_string += ("," + data[i])
					test_name = data[i-1] + '__' + data[0]
					gene_data = genes_out[genes_out_index.index(test_name)]
					for j in range(2, len(gene_data)):
						output_string += ("," + gene_data[j])						
				output_string += "\n"

		if replicons_out == [] and genes_out != []:
			data_out = []
			isolate_list = []
			for index in genes_out_index:
				name = index.split('__')
				if name[0] not in replicon_list:
					replicon_list.append(name[0])
				if name[1] not in isolate_list:
					isolate_list.append(name[1])
			if replicons != []:				
				if len(replicon_list) > 1:
					ordered_replicon_list = []
					for replicon in replicons:
						if replicon in replicon_list:
							ordered_replicon_list.append(replicon)
					replicon_list = ordered_replicon_list
			for replicon in replicon_list:
				for report in genes_list:
					if report[0] == replicon:
						for j in range(1,len(report)):
							header += (',' + report[j])
			for isolate in isolate_list:
				output_string += isolate
				for replicon in replicon_list:
					test_name = replicon + '__' + isolate
					data = genes_out[genes_out_index.index(test_name)]
					for i in range(2, len(data)):
						output_string += ("," + data[i])
				output_string += "\n"

		header += "\n"
		output_file = open(output, "w")
		output_file.write(header)		
		output_file.write(output_string)
		output_file.close()
		return
		
	### MAIN PROCESS
			
	# set up variables
	if options.input_replicon == '' and options.input_gene == '':
		print '\nNo coverage file(s) provided (-i or -I)'
		sys.exit()
	else:
		input_replicon = options.input_replicon
		input_gene = options.input_gene

	if options.output == '':
		print '\nOutput file not provided (-o)'
		sys.exit()
	else:
		output = options.output

	include_list = options.include_list

	if options.replicons != "":
		replicons = options.replicons.split(',')
	else:
		replicons = []

	try:
		test_depth = float(options.depth)
	except:
		print '\nCould not read in depth: check option -d'
		sys.exit()

	if options.coverage != "":
		try:
			test_coverage = float(options.coverage)
			if test_coverage < 0.0 or test_coverage > 100.0:
				print '\nCoverage cutoff outside expected range (0 to 100): check option -c'
				sys.exit()
		except:
			print '\nCould not read in coverage cutoff: check option -c'
			sys.exit()
	else:
		test_coverage = 0.0

	replicons_out = []
	genes_out = []
	print "Reading input file(s)..."
	replicons_out, genes_out, genes_out_index, genes_list = getHeatMapData(replicons_out, genes_out, replicons, input_replicon, input_gene, include_list, test_depth, test_coverage)
	print "\n...done\n\nWriting output file..."
	heatMapDataOut(replicons_out, genes_out, genes_out_index, genes_list, replicons, output)
	print "\n...done\n"
#	print ('Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
