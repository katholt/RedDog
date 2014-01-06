'''
dervieAllStats.py

collates the general statistics for set of reads that have gone through
the pipeline - pangenome run

example:
python deriveAllStats.py <isolate>_rep_cover.txt

Created:	29042013
Modified:	23072013
author: David Edwards
'''
import sys
from pipe_utils import splitPath

output_allStats = ""


repCoverFile = sys.argv[1]
repCover = open(repCoverFile, "r") 
(prefix, middle, ext) = splitPath(repCoverFile)
name = middle[:-10]

output_allStats += name +"\t"
name = prefix + "/" + name

#get % cover and depth for each replicon
replicon_names = []
depth_test_value = 0.0
cover = ""
depth = ""
for line in repCover:
	entry = line.split()
	replicon_names.append(entry[0])
	depth += entry[2] +"\t"
	cover += entry[3] +"\t"
output_allStats += cover + depth

#get % mapped stats for each replicon and the total
replicons_mapped = []
repMapped = open(name + "_samStats.txt")
percent_A = percent_G = percent_T = percent_C = percent_N = '0'
insert_mean = ""
insert_stdev = ""
for line in repMapped:
	if line.startswith('reads'):
		entry = line.split()
		total_reads = entry[1]
	elif line.startswith('mapped reads'):
		entry = line.split()
		mapped_reads = entry[2]
	elif line.startswith('len max'):
		entry = line.split()
		length_max = entry[2]
	elif line.startswith('insert mean'):
		entry = line.split()
		insert_mean = entry[2]
	elif line.startswith('insert stdev'):
		entry = line.split()
		insert_stdev = entry[2]
	elif line.startswith('base qual mean'):
		entry = line.split()
		base_qual_mean = entry[3]
	elif line.startswith('base qual stdev'):
		entry = line.split()
		base_qual_stdev = entry[3]
	elif line.startswith('%'):
		line = line[1:]
		if line.startswith('A\t'):
			entry = line.split()
			percent_A = entry[1]
		elif line.startswith('G\t'):
			entry = line.split()
			percent_G = entry[1]
		elif line.startswith('C\t'):
			entry = line.split()
			percent_C = entry[1]
		elif line.startswith('T\t'):
			entry = line.split()
			percent_T = entry[1]
		elif line.startswith('N\t'):
			entry =line.split()
			percent_N = entry[1]
		else:
			for replicon in replicon_names:
				if line.startswith(replicon):
					entry = line.split()
					replicons_mapped.append([replicon, entry[1]])
if insert_mean == "":
	insert_mean = "-"
if insert_stdev == "":
	insert_stdev = "-"
mapped_core_rep = 0.0
for replicon in replicon_names:
	value = ""
	for mapped in replicons_mapped:
		if replicon == mapped[0]:
			value = mapped[1]
	if value == "":
		output_allStats += '0.0\t'
	else:
		output_allStats += str((float(value)/100)*(float(mapped_reads)*100/float(total_reads))) +"\t"

output_allStats += str(float(mapped_reads)*100.0/float(total_reads)) +"\t"+ total_reads +"\t"
output_allStats += insert_mean +"\t"+ insert_stdev +"\t"+ length_max +"\t"+ base_qual_mean +"\t"+ base_qual_stdev +"\t"
output_allStats += percent_A +"\t"+ percent_T +"\t"+ percent_C +"\t"+ percent_G +"\t"+ percent_N + "\n"

out_allStats = open(name + "_AllStats.txt", "w")
out_allStats.write(output_allStats)
out_allStats.close()
