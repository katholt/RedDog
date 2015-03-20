'''
Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)
'''
# read in coverage output file (% coverage for each gene in mapped pan genome)
# remove genes that are not covered (to set % level, default 95%) in any strains
# generate presence/absence table (1/0, based on set % level, default 95%)
# generate summary of genes, reporting total number of strains with the gene
# if strain IDs and groups provided, will summarize presence by group
# if depth file is provided, will use this along with % coverage to call presence/absence (min depth = 5 by default)
# can also accept genbank annotation file and output product identifiers etc along with genes in gene summary table [adds significant time]
# can use R & rpy to plot pan and core rarefaction curve data using specified number of sampling iterations [adds significant time]
# FUTURE:
# could extract DNA/AA sequences and print these into gene summary file [increases file size significantly]
# store data in python object suitable for direct use by R, so not reliant on reading CSV data into R which is slow
# do calculations directly in python and just pass to R for plotting
import string, re
import os, sys, subprocess
import StringIO
import collections
from optparse import OptionParser

def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-g", "--gcfile", action="store", dest="gcfile", help="gene content file (csv; required)", default="")
	parser.add_option("-d", "--depthfile", action="store", dest="depthfile", help="gene depth file (csv; optional)", default="")
	parser.add_option("-i", "--info", action="store", dest="info", help="info file (optional; csv with header: id,name,group)", default="")
	parser.add_option("-r", "--ref", action="store", dest="ref", help="reference annotation (gbk; optional)", default="")
	parser.add_option("-u", "--uniqueid", action="store", dest="uniqueid", help="gene id qualifer in reference annotation (default locus_tag)", default="locus_tag")
	parser.add_option("-c", "--cutoff", action="store", dest="cutoff", help="cutoff coverage level (%) (95)", default="95")
	parser.add_option("-D", "--mindepth", action="store", dest="mindepth", help="minimum depth (5)", default="5")
	parser.add_option("-p", "--plots", action="store", dest="plots", help="generate pan-genome plots using this many iterations (default 0, ie off; module R-gcc/2.7.2 must be loaded)", default="0")
	parser.add_option("-o", "--out", action="store", dest="out", help="file to output presence/absence matrix (csv; not output if not provided)", default="")
	parser.add_option("-s", "--summary", action="store", dest="summary", help="output for genewise summary (csv; not output if not provided)", default="")
	
	return parser.parse_args()


if __name__ == "__main__":

	(options, args) = main()
	
	def checkPresence(cov,cutoff,depth,mindepth,depthfile):
		present = False # assume absent
		if (float(cov) >= float(cutoff)):
			if depthfile == "":
				present = True # no depth provided, coverage is enough
			else:
				if float(depth) > float(mindepth):
					present = True # require depth to determine presence
		return(present)
			
	header = []
	if options.out !="":
		o = file(options.out,"w")
	if options.summary !="":
		s = file(options.summary,"w")
	if options.summary =="" and options.out=="":
		DoError("No output requested. Use -o to generate presence/absence matrix and/or -s to generate gene summary file.")
				
	# read in strain info table if provided
	info_header=[]
	groups = {} # key = group, value = list of ids
	id2group = {} # key = id, value = group
	id2name = {} # key = id, value = strain name
	if options.info != "":
		f = file(options.info,"r")
		for line in f:
			if len(info_header)==0:
				info_header = line.rstrip().split(",")
			if len(info_header)>0:
				fields = line.rstrip().split(",")
				id = fields[0]
				strain = fields[1]
				id2name[id] = strain
				if len(fields) > 2:
					group = fields[2]
					id2group[id] = group
					if group in groups:
						groups[group].append(id) # add id to group
					else:
						groups[group] = [id] # new group
		f.close()
		
	# read in depth table if provided
	depth_header = []
	depth = collections.defaultdict(dict) # key 1 = gene, key 2 = strain, value = depth (float)
	if options.depthfile != "":
		f = file(options.depthfile,"r")
		for line in f:
			if len(depth_header)==0:
				depth_header = line.rstrip().split(",")
			else:
				fields = line.rstrip().split(",")
				gene = fields[0]
				for i in (range(1,len(fields))):
					id = depth_header[i]
					depth[gene][id] = float(fields[i])
		f.close()
					
	# read in genbank file, if provided
	gene_info = {} # key = gene locus tag, value = product
	if options.ref != "":
		from Bio import SeqIO
		from Bio.SeqFeature import SeqFeature, FeatureLocation
		from Bio.SeqRecord import SeqRecord
		from Bio.Alphabet import generic_dna
		for r in SeqIO.parse(options.ref, "genbank"):
			for f in r.features:
				if f.type == "CDS":
					tag = f.qualifiers[options.uniqueid][0]
					product = re.sub(',','_',f.qualifiers["product"][0]) # replace commas to prevent interfering with CSV format
					gene_info[tag] = product		
		
	numgroups = 0 # number of groups
	if len(info_header)>2:
		numgroups = len(groups) # number of groups
		
	# write header for summary table, include group names if provided
	if options.summary != "":
		s.write("gene")
		if options.ref != "":
			s.write(",product")
		s.write(",strainCount")
		if numgroups > 0:
			for group in groups:
				s.write(","+group)
		s.write("\n")
		
	f = file(options.gcfile,"r")
	gene_count = 0
	for line in f:
		if len(header)==0:
			header = line.rstrip().split(",")
			# check we have depths for these if we are using depths
			if options.depthfile != "":
				for i in (range(1,len(header))):
					if header[i] not in depth_header:
						sys.exit("ID " + header[i] + " in coverage table " + options.gcfile + " is missing from depth file " + options.depthfile)
			if options.info == "" and options.out != "":
				o.write(line) # can just write ids
			else:
				# change ids to strain names
				header_line = header[0] + "," # column 1 name
				for i in (range(1,len(header))):
					id = header[i]
					if id in id2name:
						header_line += id2name[id] + ","
					else:
						header_line += id + ","
				if options.out != "":
					o.write(header_line + "\n")
		else:
			fields = line.rstrip().split(",")
			gene = fields[0]
			# check we have depths for this gene if we are using depths
			if gene not in depth:
				DoError("gene " + gene + " in coverage table " + options.gcfile + " is missing from depth file " + options.depthfile)
			count = 0
			group_counts = {} # key = group, value = count
			for i in range(1,len(fields)):
				id = header[i]
				if checkPresence(fields[i],options.cutoff,depth[gene][id],options.mindepth,options.depthfile):
					# gene is present in this strain
					count += 1
					fields[i] = "1"
					if numgroups > 0:
						if id in id2group:
							group = id2group[id]
							if group in group_counts:
								group_counts[group] += 1
							else:
								group_counts[group] = 1
						else:
							if "ungrouped" in group_counts:
								group_counts["ungrouped"] += 1
							else:
								group_counts["ungrouped"] = 1	
				else:
					fields[i] = "0" # gene is not present in this strain
			if (count > 0):
				# otherwise exclude gene from outputs
				gene_count += 1
				if options.out !="":
					o.write((",").join(fields)+"\n")
				if options.summary !="":
					s.write(gene)
					if options.ref != "":
						s.write(","+gene_info[gene])
					s.write(","+str(count))
					if numgroups>0:
						for group in groups:
							if group in group_counts:
								s.write(","+str(float(group_counts[group])/float(len(groups[group]))))
							else:
								s.write(",0") # none encountered in this group
					s.write("\n")
	if options.out != "":
		o.close()
	if options.summary != "":
		s.close()
	
	# generate pan genome plots using Rpy
	ylim = str((gene_count / 1000 + 1) * 1000)
	if options.plots != "0" and options.out != "":
		from rpy import *	
		r('content<-read.csv("'+options.out+'",row.names=1)')
		r('content_bygene<-table(apply(content,1,sum))')
		r('content0<-0')
		r('if (0 %in% names(content_bygene)) {content0<-content_bygene[names(content_bygene)==0]}')
		r('content0N<-0')
		r('if (ncol(content) %in% names(content_bygene)) {content0N<-content_bygene[names(content_bygene)==ncol(content)]}')
		r('counts_core<-matrix(ncol=2)')
		r('counts_core[1,]<-c(ncol(content),content0N)')
		r('counts_pan<-matrix(ncol=2)')
		r('counts_pan[1,]<-c(ncol(content),nrow(content)-content0)')		
		r('for (i in 2:(ncol(content)-1)) {for (j in 1:' + options.plots + ') {content_sub<-content[,round(runif(i,1,ncol(content)))] ; content_sum<-apply(content_sub,1,sum);content_core<-length(content_sum[content_sum==i]);content_pan<-length(content_sum[content_sum>0]);counts_core<-rbind(counts_core,c(i,content_core));counts_pan<-rbind(counts_pan,c(i,content_pan))}}')
		r('write.table(counts_core,file="' + options.out + '_' + options.plots + '_CoreGenes.csv",sep=",")')
		r('write.table(counts_pan,file="' + options.out + '_' + options.plots + '_PanGenome.csv",sep=",")')
		r('csh<-boxplot(counts_core[,2]~counts_core[,1])')
		r('psh<-boxplot(counts_pan[,2]~counts_pan[,1])')
		r('pdf(file="' + options.out + '_' + options.plots + '_PanGenome.pdf")')
		r('plot(psh$names,psh$stats[3,],ylim=c(0,'+ylim+'),pch="",main="Pan genome",xlab="N isolates",ylab="N genes")')
		r('lines(psh$names,psh$stats[3,],lwd=2)')
		r('lines(psh$names,psh$stats[2,],lty=2)')
		r('lines(psh$names,psh$stats[4,],lty=2)')
		r('lines(csh$names,csh$stats[3,],lwd=2)')
		r('lines(csh$names,csh$stats[2,],lty=2)')
		r('lines(csh$names,csh$stats[4,],lty=2)')
		r('dev.off()')
				