'''
Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)
'''
#!/usr/bin/env python
#
# Read in SNP alleles (csv format)
#  can take a file containing a list of strains to include (-l), otherwise all are included
#  can take a list of outgroups (-o) => if specified SNPs that are non-variable in in the ingroup are removed on read-in; and variation in the outgroups is ignored in assessing conservation
#  Note it doesn't matter whether or not the outgroup is included in the list of strains.
# The table is then parsed as specified by the modules in -m, these include:
#   aln - convert to fasta alignment
#   filter - filter SNPs that are included/excluded in regions specified via -x (genbank, gff or 2-column CSV table format)
#   clean - filter out any pairs of SNPs with -P bp between them (default 3bp, minimum 2bp) and any trio or more of SNPs within -W bp in any isolate (default 10bp, minimum 3bp or -P if greater than 3)
#   cons - filter SNP positions that are not conserved above a cutoff specified via -c (e.g. -c 0.99 -> all snps with >1% missing alleles is filtered out)
#   core - filter SNPs not in genes that are conserved with % coverage cutoff, specified via -Z, obtained from the gene cover table, specified via -z, across all core isolates, specified via -L. As for -l, outgroups are ignored.   
#   vcf - convert snp table to vcf format (requires -v name of reference mapped, optional -V name of isolate mapped, otherwise -v used)
#         can also produce VCF with SNPs filtered out at each stage (set '-A True' and add vcf after each filtering step)
#   fasttree - submit a fasttree job to SLURM
#   rax - submit threaded RAxML jobs to SLURM (can specify walltime, memory, threads, number of boostraps and number of replicate runs)
# Any number of these modules can be supplied in any order; the order they are given is the order they will be run
#   specify modules in a comma-separated list, e.g. '-m filter,cons,aln,rax' will run region filtering, then conservation filter, then make a fasta alignment and use this as input to 5 raxml jobs with 100 bootstraps using 8 threads
#   clean should only be run after any recombinant regions have been indentified for exclusion
#
# This version also handles multiple sequence genbank files for coding entries. The sequence in the genbank relevant to the SNP table must be specified (-q queryseq).  
# It is also quicker generating coding consequences, and now uses much less memory (about twice the size of the allele table) 

# NOTE this script submits fasttree or RAxML jobs to SLURM
#
# Authors - Kat holt (kholt@unimelb.edu.au)
#         - David Edwards (d.edwards2@student.unimelb.edu.au) 
#
# Example command on barcoo:
'''
module load python-gcc/2.7.5
python /vlsci/VR0082/shared/code/holtlab/parseSNPtable_multigbk.py -s snps.csv -p prefix -r genbank -q queryseq -m aln,coding,rax
'''
#
# Last modified - Oct 7, 2014
# Changes:
#	 15/10/13 - added strain subset option
#    25/03/14 - added multiple sequence genbank file handling 
#             - improved 'coding' option performance
#    27/05/14 - changes to improve memory performance and filter of regions (especially overlapping regions)
#			  - also added -d directory option
#		08/14 - update to RAXML submission
#	 12/09/14 - fix for filter of regions: filtered table now passed back correctly
#			  - also updated inital SNP table reading message
#	 28/09/14 - added cleaning (filtering) of erroneous SNPs
#    01/10/14 - added output of SNP table to vcf format
#    03/10/14 - fixed minor error in output of compound vcf format for Gingr
#    07/10/14 - added filtering for core SNPs as specified by Gene Coverage table from RedDog

import os, sys, subprocess, string, re, random
import collections
import operator
from optparse import OptionParser
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
#import resource

def main():
	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	# output options
	parser.add_option("-p", "--prefix", action="store", dest="prefix", help="prefix to add to output files (default none)", default="")
	parser.add_option("-d", "--directory", action="store", dest="directory", help="directory to send output files (default none)", default="")
	
	# modules to run
	parser.add_option("-m", "--modules", action="store", dest="modules", help="modules to run, comma separated list in order (default aln, e.g. filter,vcf,cons,vcf,aln,fasttree,rax,core,vcf,aln,fasttree,rax,coding)", default="aln")

	# snptable reading
	parser.add_option("-s", "--snptable", action="store", dest="snptable", help="SNP table (CSV)", default="")
	parser.add_option("-g", "--gapchar", action="store", dest="gapchar", help="gap character (default -, could be N)", default="-")
	parser.add_option("-o", "--outgroup", action="store", dest="outgroup", help="comma separated list; outgroup strains (alleles will be included but not sites that vary only in outgroups)", default="")
	parser.add_option("-l", "--subset", action="store", dest="subset", help="file containing list of strains to include (one per line), otherwise all strains included", default="")

	# region filtering
	parser.add_option("-x", "--regions", action="store", dest="regions", help="file of regions to include/exclude (gbk)", default="")
	parser.add_option("-y", "--include", action="store", dest="include", help="include (default exclude)", default="exclude")

	# conservation filtering
	parser.add_option("-c", "--conservation", action="store", dest="conservation", help="minimum conservation across samples required to retain SNP locus (default 0.99)", default="0.99")

	# snp cleaning
	parser.add_option("-P", "--pairs", action="store", dest="pairs", help="maximum distance between pairs of SNPs to remove (default 3bp, minimum 2bp)", default="3")
	parser.add_option("-W", "--window", action="store", dest="window", help="snp window to check clusters with more than three SNPs (default 10bp, minimum 3bp or -P if greater than 3)", default="10")	
	
	# core gene filtering and/or coding consequences
	parser.add_option("-r", "--refseq", action="store", dest="refseq", help="reference sequence file (gbk)", default="")
	parser.add_option("-q", "--queryseq", action="store", dest="queryseq", help="query sequence in reference sequence file (multisequence gbk)", default="")

	# core gene filtering
	parser.add_option("-L", "--core_strains", action="store", dest="core_strains", help="file containing list of strains to include in core genome (one per line, outgroup[s] ignored), otherwise all strains sans outgroup(s) included", default="")
	parser.add_option("-z", "--gene_coverage", action="store", dest="gene_coverage", help="gene coverage table (CSV)", default="")
	parser.add_option("-Z", "--core_coverage", action="store", dest="core_coverage", help="minimum % coverage of each gene (as ratio) across core_isolates required to retain SNP locus (default 0.90 - 90%)", default="0.9")

	# coding consequences
	parser.add_option("-f", "--genefeatures", action="store", dest="genefeatures", help="feature types for protein coding genes (default CDS; can be multiple comma-sep)", default="CDS")
	parser.add_option("-e", "--excludefeatures", action="store", dest="excludefeatures", help="feature types to exclude (default gene,misc_feature)", default="gene,misc_feature")
	parser.add_option("-i", "--identifier", action="store", dest="identifier", help="unique identifier for features (locus_tag)", default="locus_tag")

    # snptable to vcf
	parser.add_option("-v", "--reference_name", action="store", dest="reference_name", help="name of reference mapped (required, use MT for PLINK", default="")
	parser.add_option("-V", "--ref_isolate_name", action="store", dest="ref_isolate_name", help="name of the isolate used for reference (default reference_name)", default="")
	parser.add_option("-A", "--add_vcf", action="store_true", dest="add_vcf", help="Set to get VCF with PASS/FAIL for all SNPs in first SNP table [post-isolate/outgroup filtering] (default False, need at least one 'vcf' module assigned to -m)", default=False)

	# rax parameters
	parser.add_option("-t", "--walltime", action="store", dest="walltime", help="walltime for raxml jobs (default 5-12:0, ie 5.5. days)", default="5-12:0")
	parser.add_option("-M", "--memory", action="store", dest="memory", help="memory for raxml jobs (24576)", default="24576")
	parser.add_option("-T", "--threads", action="store", dest="threads", help="number of threads per raxml job (8)", default="8")
	parser.add_option("-N", "--N", action="store", dest="N", help="number of raxml bootstraps (100)", default="100")
	parser.add_option("-n", "--numrax", action="store", dest="numrax", help="number of raxml jobs (5)", default="5")

	return parser.parse_args()

if __name__ == "__main__":

	(options, args) = main()
	nt = ["A","C","G","T"]

	def isVar(alleles):
		numAlleles = 0
		for a in nt:
			if alleles.count(a) > 0:
				numAlleles += 1
		return (numAlleles > 1)
	
	# read csv; return as dictionary of dictionaries and list of strains
	def readSNPTable(infile,outgroup_list,strain_list_file,pre):
		print "\nReading SNP table from " + infile
		outgroups = [] # list of outgroups provided
		outgroups_used = [] # list of outgroups encountered
		if options.outgroup != "":
			outgroups = options.outgroup.split(",")
			print "    outgroup(s): " + ",".join(outgroups)
			pre += "_" + str(len(outgroups)) + "outgroup"
			if len(outgroups) > 1:
				pre += "s"
				
		strainlist = [] # list of strains to include, excluding outgroups
		if strain_list_file != "":
			f = file(strain_list_file,"r")
			for line in f:
				strain = line.rstrip()
				if strain not in outgroups:
					strainlist.append(strain)
			f.close()
			print "    including " + str(len(strainlist)) + " ingroup strains listed in file " + strain_list_file
			pre += "_" + str(len(strainlist)) + "strains"
		else:
			print "    including all strains"
		snptable = []
		strains = [] # strains from header
		ignored = []
		pre += "_var"
		f = file(infile, "r")
		lines = f.readlines()
		f.close()
		o = file(pre + ".csv","w")
		o.write("Pos")
		fields = []
		count = 0
		for i in range(len(lines)):
#			if i % 10000 == 0:
#				print 'Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
			if fields == []:
				fields = lines[i].rstrip().split(',')
			if len(strains)==0:
				strains = fields
				if len(strainlist) == 0:
					for j in range(1,len(strains)):
						if strains[j] not in outgroups:
							strainlist.append(strains[j]) # retain all strains (except outgroups)
				# remove strains from the strainlist if we have not encountered them in the actual table
				for strain in strainlist:
					if strain not in strains:
						strainlist.remove(strain)
				# print header for new table
				for j in range(1, len(fields)):
					if strains[j] in strainlist or strains[j] in outgroups:
						o.write(","+strains[j])						
				o.write("\n")
			else:
				j=0
				snp = ''
				while lines[i][j] != ',':
					snp += lines[i][j]
					j+=1
				alleles_ingroup = [] # alleles for this snp
				alleles_outgroup = [] # alleles for outgroups
				for k in range(1,len(fields)):
					if strains[k] in strainlist and strains[k] not in outgroups and lines[i][(j+2*(k-1)+1)] not in alleles_ingroup:
						alleles_ingroup.append(lines[i][(j+2*(k-1)+1)]) # add this allele to the list
					elif strains[k] in outgroups and lines[i][(j+2*(k-1)+1)] not in alleles_outgroup:
						alleles_outgroup.append(lines[i][(j+2*(k-1)+1)]) # add this allele to the outgroup list
						if strains[k] not in outgroups_used:
							outgroups_used.append(strains[k])
				if isVar(alleles_ingroup):
					o.write(str(snp))
					snptable.append([])
					snptable[count].append(snp)
					snptable[count].append('')
					# variable in ingroup, store alleles for included strains
					for k in range(1, len(fields)):
						if strains[k] in strainlist or strains[k] in outgroups:
							snptable[count][1] += lines[i][(j+2*(k-1)+1)].upper()
							o.write(","+lines[i][(j+2*(k-1)+1)].upper())
					o.write("\n")
					count +=1
				else:
					ignored.append(snp)
		o.close()
		strains.pop(0) # remove SNP column header
		strains_used = []
		for strain in strains:
			if strain in strainlist or strain in outgroups_used:
				strains_used.append(strain) 			

		print "\n... finished reading " + str(len(snptable) + len(ignored)) + " SNPs in total"
		print "... keeping " + str(len(snptable)) + " variable SNPs in " + str(len(strainlist)) + " ingroup strains"
		print "... ignoring " + str(len(ignored)) + " SNPs that are non-variable among these ingroup strains"
		return(snptable, strains_used, pre) # include outgroups that appear in the snptable in strainlist
	
	def printFasta(snptable, strains, outfile):
		print "\nPrinting alignment to file " + outfile
		o = file(outfile,"w")
		for strain in range(len(strains)):
			o.write(">" + strains[strain] + "\n")
			seq = ""
			for snp in range(len(snptable)): # cycle over SNPs
				o.write(str(snptable[snp][1][strain]))
			o.write("\n")
		o.close()
		print "\n... done"

	def getCodons(genestart,genestop,genestrand,snp,derived,ancestral,sequence,complement):
		codon = ()
		posincodon = 0
		# determine coordinates of codon within genome
		if genestrand == 1:
			posingene = snp-genestart # note genestart is in -1 offset space, snp is not
			if posingene % 3 == 0:
				codon = (snp-2,snp-1,snp)
				posincodon = 3
			elif posingene % 3 == 1:
				codon = (snp,snp+1,snp+2)
				posincodon = 1
			else:
				codon = (snp-1,snp,snp+1)
				posincodon = 2
		elif genestrand == -1:
			posingene = genestop-snp+1 # note genestop is not in -1 offset space
			if posingene % 3 == 0:
				codon = (snp+2,snp+1,snp)
				posincodon = 3
			elif posingene % 3 == 1:
				codon = (snp,snp-1,snp-2)
				posincodon = 1
			else:
				codon = (snp+1,snp,snp-1)
				posincodon = 2
		else:
			DoError("Unrecognised gene strand:" + genestrand)
		# extract codon sequence from reference genome
		codonseq = [ str(sequence[codon[0]-1]), str(sequence[codon[1]-1]) , str(sequence[codon[2]-1]) ] # codon sequence
		if genestrand == -1:
			# complement the reverse strand
			codonseq = [ complement[codonseq[0]], complement[codonseq[1]] , complement[codonseq[2]] ] # codon sequence
		# insert ancestral base
		if genestrand == 1:
			codonseq[posincodon-1] = ancestral # replace snp within codon
		elif genestrand == -1:
			codonseq[posincodon-1] = complement[ancestral] # replace snp within codon
		ancestral_codon = Seq(''.join(codonseq),IUPAC.unambiguous_dna)
		# mutate with current SNP
		if genestrand == 1:
			codonseq[posincodon-1] = derived # replace snp within codon
		elif genestrand == -1:
			codonseq[posincodon-1] = complement[derived] # replace snp within codon
		derived_codon = Seq(''.join(codonseq),IUPAC.unambiguous_dna)
		ancestralAA = ancestral_codon.translate()
		derivedAA = derived_codon.translate()
		return(ancestral_codon,derived_codon,ancestralAA,derivedAA,posingene,posincodon)

	def runCoding(pre, snptable, options, complement):
		if options.refseq=="":
			print "\nNo reference genbank file specified (-r), can't do coding analysis"
		else:
			genefeatures = options.genefeatures.split(",")
			excludefeatures = options.excludefeatures.split(",")
			# order SNPs
			snp_list_ordered = []
			snp_list_paired = []
			for snp in range(len(snptable)): # cycle over SNPs
				snp_list_paired.append([snp,int(snptable[snp][0])])
			snp_list_paired.sort(key=operator.itemgetter(1))
			for snp in range(len(snptable)):
				snp_list_ordered.append(snp_list_paired[snp][0])
			print "\nReading gene features from reference " + options.refseq
			# check coding consequences and generate genbank file of SNP loci
			## READ IN GENBANK FILE 
			Passed = True
			handle = open(options.refseq,"r")
			if options.queryseq=="":
				try:
					record = SeqIO.read(handle, "genbank")
					sequence = record.seq
					geneannot = record.features
				except:
					Passed = False
					print "\nCheck reference sequence for multiple records: can't do coding analysis"
			else:
				records = SeqIO.parse(handle, "genbank")
				Passed = False
				for item in records:
					if item.name==options.queryseq:
						record = item
						sequence = SeqRecord(item.seq)
						geneannot = item.features
						Passed = True
				if Passed == False:
					print "\nCheck reference sequence: queryseq (-q) not found"
			if Passed==True:
				print "Determining coding changes"
				## GET CONSEQUENCES FOR SNPS and WRITE SNP ANNOTATION FILE
				# first make index for features
				feature_list = []
				feature_count = 0
				for feature in geneannot:
					if feature.type != "source" and feature.type not in excludefeatures:
						start = feature.location.nofuzzy_start
						stop = feature.location.nofuzzy_end
						feature_list.append([start,stop,feature_count])
					feature_count += 1
				feature_slice = []
				if len(feature_list) > 0:
					slice_size = len(sequence)/len(feature_list)+1
					for slice in range((len(sequence)/slice_size)+2):
						feature_slice.append([])
				else:
					slice_size = len(record) +1
					feature_slice.append([])
					feature_slice.append([])
				feature_count=0
				for feature in feature_list:
					slice1 = feature_list[feature_count][0]/slice_size
					slice2 = feature_list[feature_count][1]/slice_size
					feature_slice[slice1].append([feature_list[feature_count][0],feature_list[feature_count][1],feature_list[feature_count][2]])
					while slice1 < slice2:
						slice1 += 1
						feature_slice[slice1].append([feature_list[feature_count][0],feature_list[feature_count][1],feature_list[feature_count][2]])
					feature_count += 1
				o = open(pre + "_consequences.txt","w")
				o.write("\t".join(["SNP","ref","alt","change","gene","ancestralCodon","derivedCodon","ancestralAA","derivedAA","product","ntInGene","codonInGene","posInCodon","\n"])) # header
				intergenic_count = 0
				ns_count = 0
				syn_count = 0
				other_feature_count = 0
				for snp in snp_list_ordered:
					ref_allele = sequence[int(snptable[snp][0])-1]
					allele_list = []
					for strain in range(len(snptable[snp][1])):
						if snptable[snp][1][strain] not in allele_list:
							allele_list.append(snptable[snp][1][strain])
					if options.gapchar in allele_list:
						allele_list.remove(options.gapchar)
					if ref_allele in allele_list:
						allele_list.remove(ref_allele)
					if len(allele_list)>0:
						for alt_allele in allele_list:
							hit = 0 # initialize
							snp_slice = int(snptable[snp][0])/slice_size
							if feature_slice[snp_slice] != []:
								for feature_index in feature_slice[snp_slice]:
									if int(snptable[snp][0]) in geneannot[feature_index[2]]:
										hit = 1
										start = int(geneannot[feature_index[2]].location.nofuzzy_start) # feature start
										stop = int(geneannot[feature_index[2]].location.nofuzzy_end) # feature stop
										id = ""
										product = ""
										if options.identifier in geneannot[feature_index[2]].qualifiers:
											id = geneannot[feature_index[2]].qualifiers[options.identifier][0]
										if 'product' in geneannot[feature_index[2]].qualifiers:
											product = geneannot[feature_index[2]].qualifiers['product'][0]
										if geneannot[feature_index[2]].type in genefeatures:
											# get coding effect of coding features
											(ancestral_codon,derived_codon,ancestralAA,derivedAA,posingene,posincodon)=getCodons(start,stop,geneannot[feature_index[2]].strand,int(snptable[snp][0]),alt_allele,ref_allele,sequence,complement)
											change = "ns"
											if str(ancestralAA) == str(derivedAA):
												change = "s"
												syn_count += 1
											else:
												ns_count += 1
											# add SNP to genbank
											codon_number = posingene / 3
											if posincodon != 3:
												codon_number += 1
											note = change + " SNP " + ref_allele + "->" + alt_allele + " at nt " + str(posingene) + ", position " + str(posincodon) + " in codon " + str(codon_number) + "; " + str(ancestral_codon) + "->" + str(derived_codon) + "; " + str(ancestralAA) + "->" + str(derivedAA) 
											record.features.append(SeqFeature(FeatureLocation(int(snptable[snp][0])-1,int(snptable[snp][0])), type="variation", strand=1, qualifiers = {'note' : [note]}))
											o.write("\t".join([snptable[snp][0],ref_allele,alt_allele,change,id,str(ancestral_codon),str(derived_codon),str(ancestralAA),str(derivedAA),product,str(posingene),str(codon_number),str(posincodon),"\n"]))
										else:
											# non-protein coding feature
											other_feature_count += 1
											if geneannot[feature_index[2]].strand == 1:
												posingene = int(snptable[snp][0])-start+1
											else:
												posingene = stop-int(snptable[snp][0])+1
											o.write("\t".join([snptable[snp][0],str(ref_allele),str(alt_allele),geneannot[feature_index[2]].type,id,"","","","",product,str(posingene),"","","\n"]))
											record.features.append(SeqFeature(FeatureLocation(int(snptable[snp][0])-1,int(snptable[snp][0])), type="variation", strand=1, qualifiers = {'note' : ["SNP " + ref_allele + "->" + alt_allele + " in non-CDS feature" ]}))
							if hit == 0:
								# SNP is intergenic
								intergenic_count += 1
								o.write("\t".join([snptable[snp][0],str(ref_allele),str(alt_allele),"intergenic","","","","","","","","\n"]))
								record.features.append(SeqFeature(FeatureLocation(int(snptable[snp][0])-1,int(snptable[snp][0])), type="variation", strand=1, qualifiers = {'note' : ["intergenic SNP " + ref_allele + "->" + alt_allele]}))
				o.close()
				SeqIO.write(record, pre + ".gbk", "genbank")
				print "\n... " + str(ns_count) + " nonsyonymous, " + str(syn_count) + " synonymous, " + str(other_feature_count) + " in other features, " + str(intergenic_count) + " in non-coding regions"
				print "... coding consequences written to file " + pre + "_consequences.txt"
				print "... SNP loci annotated in genbank file " + pre + ".gbk"

	def runFasttree(pre, snptable, strains):
		aln = pre + ".mfasta"
		if not os.path.exists(aln):
			printFasta(snptable, strains, aln) # make alignment first
		jobscript = pre + "_FastTree.sh"
		o = file(jobscript, "w")
		print "\nRunning fasttree on " + aln + ", using job script: " + jobscript
		o.write("#!/bin/bash\n#SBATCH -p main\n##SBATCH --job-name='ft" + pre)
		o.write("'\n#SBATCH --time=0-24:0")
		o.write("\n#SBATCH --mem-per-cpu=24576")
		o.write("\n#SBATCH --ntasks=1")
		o.write("\ncd " + os.getcwd())
		o.write("\nmodule load fasttree-intel\n")
		o.write("FastTree -gtr -gamma -nt " + aln + " > " + pre + ".tree\n")
		o.close()
		os.system('sbatch ' + jobscript)
		print "\n... output tree will be in " + pre + ".tree"
		
	def runRax(pre, options, snptable, strains):
		aln = pre + ".mfasta"
		if not os.path.exists(aln):
			printFasta(snptable, strains, aln) # make alignment first
		for rep in range(0,int(options.numrax)):
			# get random seeds
			seed = random.randint(10000,1000000)
			if seed % 2 == 0:
				seed += 1
			p = random.randint(10000,1000000)
			if p % 2 == 0:
				p += 1	
			# prepare job script
			jobscript = pre + "_rax_" + str(rep) + ".sh"
			o = file(jobscript, "w")
			print "\nRunning RAxML 7.7.2 (PTHREADS-SSE3) on " + aln + ", using job script: " + jobscript
			rax_pre = pre + "_" + str(rep)
			if os.path.exists("RAxML_info."+rax_pre):
				rax_pre = pre + "_" + str(rep) + "_" + str(seed) # make sure output is unique
			o.write("#!/bin/bash")
			o.write("\n#SBATCH -p main")
			o.write("\n#SBATCH --job-name=rax" + pre + str(rep))
			o.write("\n#SBATCH --time=" + options.walltime)
			o.write("\n#SBATCH --mem=" + options.memory)
			o.write("\n#SBATCH --ntasks=" + options.threads)
			o.write("\n#SBATCH --nodes=1")
			o.write("\n#SBATCH --exclusive")
			o.write("\ncd " + os.getcwd())
			o.write("\nmodule load raxml-intel/7.7.2\n")
			#o.write("raxmlHPC -s " + aln)
			o.write("raxmlHPC-PTHREADS-SSE3 -T " + options.threads + " -s " + aln)
			o.write(" -n " + rax_pre + " -f a -m GTRGAMMA -x " + str(seed) )
			o.write(" -N " + options.N + " -p " + str(p) + "\n")
			o.close()
			os.system('sbatch ' + jobscript)
			print "\n... output will be in RAxML*" + rax_pre
			
	def filter(snptable, strainlist, pre, options):	
		# parse genomic regions
		if options.regions =="":
			print "\nNo regions file provided (-x), can't filter SNPs from genomic regions"
		else:
			regions = [] # list of positions to exclude/include
			(f,ext) = os.path.splitext(options.regions)
			if ext == ".gbk" or ext == ".gb":
				handle = open(options.regions,"r")
				record = SeqIO.read(handle, "genbank")
				for region in record.features:
					region_start = min(int(region.location.nofuzzy_start),int(region.location.nofuzzy_end))
					region_stop = max(int(region.location.nofuzzy_start),int(region.location.nofuzzy_end))
					regions.append([region_start,region_stop])
			else: 
				start = 0
				stop = 1
				delim = "," # assume csv
				if ext == ".gff":
					start = 3
					stop = 4
					delim = "\t" # gff is tab-delim
				elif ext == ".txt":
					# assume text table exported from genbank
					start = 2
					stop = 3
					delim = "\t"
				e = file(options.regions, "r")
				for line in e:
					fields = line.rstrip().split(delim)
					region_start = min(int(fields[start]),int(fields[stop]))						
					region_stop = max(int(fields[start]),int(fields[stop]))
					if region_start == int(fields[stop]):
						region_start += 1
						region_stop += 1
					regions.append([region_start,region_stop])
				e.close()
			region_list = []
			regions.sort(key=operator.itemgetter(1))
			regions.sort(key=operator.itemgetter(0))
			x = 0
			while x < len(regions):
				if x == len(regions) - 1: #last region in list
					region_list.append(regions[x])
				elif regions[x][0] < regions[x+1][0] and regions[x][1] < regions[x+1][0]:
					region_list.append(regions[x])
				elif regions[x][1] >= regions[x+1][0]:
					if regions[x][1] <= regions[x+1][1]:
						regions[x+1] = [regions[x][0], regions[x+1][1]]
					else:
						regions[x+1] = regions[x]
				x += 1
			total_masked_bases = 0
			for region in region_list:
				total_masked_bases += region[1] - region[0]
			last_region_call = region_list[-1][1]
			region_slice = []
			if len(region_list) > 0:
				slice_size = (last_region_call + 1)/len(region_list)+1
				for slice in range(((last_region_call + 1)/slice_size)+1):
					region_slice.append([])
			else:
				slice_size = (last_region_call + 1) +1
				region_slice.append([])
				region_slice.append([])
			for region in region_list:
				slice1 = region[0]/slice_size
				slice2 = region[1]/slice_size
				region_slice[slice1].append([region[0],region[1]])
				while slice1 < slice2:
					slice1 += 1
					region_slice[slice1].append([region[0],region[1]])
			# filter snps
			o = file(pre + "_regionFiltered.csv","w")
			o.write(",".join(["Pos"] + strainlist) + "\n") # write header
			print "\nFiltering SNPs that are ",
			if options.include != "include":
				# excluding snps 
				print "located in excluded regions",
			else:
				print "located outside included regions",
			print "totalling " + str(total_masked_bases) + " bases\n",
			print "specified in file " + options.regions
			snpcount = 0
			to_remove = [] # list of snps to remove
			for snp in snptable:
				snp_slice = int(snp[0])/slice_size
				if options.include != "include":
					if int(snp[0]) > last_region_call:
						keep = True
					elif region_slice[snp_slice] == []:
						keep = True
					else:
						keep = True
						for region in region_slice[snp_slice]:
							if int(snp[0]) >= region[0] and int(snp[0]) < region[1]:
								keep = False
				elif options.include == "include":
					if int(snp[0]) > last_region_call:
						keep = False
					elif region_slice[snp_slice] == []:
						keep = False
					else:
						keep = False
						for region in region_slice[snp_slice]:
							if int(snp[0]) >= region[0] and int(snp[0]) < region[1]:
								keep = True
				if keep:
					o.write(snp[0])
					for strain in range(len(strainlist)):
						o.write(","+snp[1][strain])
					o.write("\n")
					snpcount += 1
				else:
					to_remove.append(int(snp[0]))
			o.close()
			pre += "_regionFiltered"
			print "... " + str(snpcount) + " SNPs passed filter; printed to " + pre + ".csv"

			if to_remove != []:
				filtered_snptable = []
				for snp in range(len(snptable)):
					if int(snptable[snp][0]) not in to_remove:
						filtered_snptable.append(snptable[snp])
				snptable = filtered_snptable
		return pre, snptable
	
	def filterCons(snptable, strainlist, pre, options, outgroups):
		try:
			cons = float(options.conservation)
			print "\nFiltering SNPs with fewer than " + str(round(100*float(options.conservation),5)) + "% known alleles",
			if len(outgroups) > 0:
				print "amongst ingroups"
			pre += "_cons" + options.conservation
			outfile = pre + ".csv"
			o = open(outfile,"w")	
			o.write(",".join(["Pos"] + strainlist) + "\n") # write header
			to_remove = [] # list of snps to remove
			snpcount = 0
			for snp in range(len(snptable)):
				n = len(strainlist) # total strains
				alleles = []
				for strain in range(len(strainlist)):
					alleles.append(snptable[snp][1][strain])
				if len(outgroups) == 0:
					allele_list = alleles # all alleles
				else:
					allele_list = []
					for strain in range(len(strainlist)):
						if strainlist[strain] not in outgroups:
							allele_list.append(alleles[strain])
						else:
							n -= 1 #correct n for outgroup
				numgaps = allele_list.count(options.gapchar) # number of gaps at this position, excluding outgroups
				callrate = 1 - float(numgaps) / n
				if callrate >= float(options.conservation):
					o.write(snptable[snp][0])
					for strain in range(len(strainlist)):
						o.write(","+alleles[strain])
					o.write("\n")
					snpcount += 1
				else:
					to_remove.append(snp) # remove SNP from table
			o.close()
			print "\n... " + str(snpcount) + " SNPs passed filter; printed to " + pre + ".csv"

			if to_remove != []:
				filtered_snptable = []
				for snp in range(len(snptable)):
					if snp not in to_remove:
						filtered_snptable.append(snptable[snp])
				snptable = filtered_snptable
		except:
			print "\nCouldn't filter SNPs based on missing alleles, couldn't understand the proportion given: -c " + options.conservation
		return pre, snptable

	def cleanSNPs(snptable, strainlist, pre, options):
		try:
			pairs = int(options.pairs)
			if pairs < 2:
				pairs = 2
			window = int(options.window)
			if window < pairs:
				window = pairs
			if window == 2:
				window += 1
			print "\nFiltering SNP pairs within " + str(pairs) + "bp (minimum 2bp)"
			print "Also filtering when 3 or more SNPs found within window of " + str(window) + "bp (minimum 3 bp or -P, if greater than 3)"
			to_remove = []
			for j in range(1,len(strainlist)):
				all_snps = [] #all snps for one strain
				for i in range(len(snptable)-1):
					if snptable[i][1][j] != snptable[i][1][0]:
						if snptable[i][1][j] in nt:
							all_snps.append(i)
				if len(all_snps) >= 2: #check pairs of snps
					for i in range(len(all_snps)-1):
						if  int(snptable[all_snps[i+1]][0]) - int(snptable[all_snps[i]][0]) < pairs:
							if all_snps[i] not in to_remove:
								to_remove.append(all_snps[i])
							if all_snps[i+1] not in to_remove:
								to_remove.append(all_snps[i+1])
				if len(all_snps) >= 3: #check triplets of SNPs (or more) within window
					for i in range(len(all_snps)-2):
						a = 1
						while (i+a < len(all_snps)) and (int(snptable[all_snps[i+a]][0]) - int(snptable[all_snps[i]][0]) < window):
							a += 1
						if a > 2:
							for b in range(a):
								if all_snps[i+b] not in to_remove:
									to_remove.append(all_snps[i+b])
			if to_remove != []:
				filtered_snptable = []
				for snp in range(len(snptable)):
					if snp not in to_remove:
						filtered_snptable.append(snptable[snp])
				snptable = filtered_snptable
			pre += "_cleanP" + str(pairs) + "W" + str(window)
			outfile = pre + ".csv"
			o = open(outfile,"w")	
			o.write(",".join(["Pos"] + strainlist) + "\n") # write header
			for snp in range(len(snptable)):
				o.write(snptable[snp][0])
				for strain in range(len(strainlist)):
					o.write(","+snptable[snp][1][strain])
				o.write("\n")
			print "\n... " + str(len(to_remove)) + " SNPs failed one or both filters"
			print "... " + str(len(snptable)) + " SNPs passed filters; printed to " + pre + ".csv"
		except:
			print "\nCouldn't filter SNPs based on pairs/window, couldn't understand the option(s) given: -P " + options.pairs + ", and/or -W " + options.window
		return pre, snptable

	def printVCF(snptable, strains, outfile, options, add_vcf):
		print "\nPrinting alignment to VCF file " + outfile + ".vcf"
		try:
			if options.reference_name == "":
				print "\nNo reference name provided, can't create VCF"
				if add_vcf:
					return ""				
			else:
				reference_name = options.reference_name
				if options.ref_isolate_name == "":
					ref_isolate_name = options.reference_name
				vcf_output = ('##fileformat=VCFv4.1\n##source=parseSNPtable '+ outfile +'\n'
					+'##INFO=<ID=WT,Number=1,Type=Integer,Description="Number of samples called reference (wild-type)">\n' 
					+'##INFO=<ID=HOM,Number=1,Type=Integer,Description="Number of samples called variant">\n'
					+'##INFO=<ID=NC,Number=1,Type=Integer,Description="Number of samples not called">\n'
					+'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
					+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+ref_isolate_name)
				outfile += ".vcf"
				for i in range(1,len(strains)):
					vcf_output += ("\t" + strains[i])
				vcf_output += "\n"
				for i in range(len(snptable)): # cycle over SNPs
					output_line = (reference_name +'\t'+str(snptable[i][0])+'\t.\t')
					variants = []
					ref_count = 0
					alt_count = 0
					null_count = 0
					for j in range(len(snptable[i][1])):
						if snptable[i][1][j] in nt:
							if snptable[i][1][j] not in variants:
								variants.append(snptable[i][1][j])
					output_line += (variants[0] + '\t')
					if len(variants) > 1:
						output_line += variants[1]
					elif len(variants) == 1:
						output_line += '.'
					if len(variants) > 2:
						for alt in range(2, len(variants)):
							output_line += (','+variants[alt])
					output_line += '\t.\tPASS\t'
					for j in range(len(snptable[i][1])):
						if snptable[i][1][j] not in variants:
							null_count += 1
						elif snptable[i][1][j] != variants[0]:
							alt_count += 1
						else:
							ref_count += 1
					output_line += ('WT='+ str(ref_count) +';HOM='+ str(alt_count) +';NC='+ str(null_count) +'\tGT')
					for j in range(len(snptable[i][1])):
						if snptable[i][1][j] not in variants:
							output_line += '\t.'
						else:
							output_line += ('\t' + str(variants.index(snptable[i][1][j])))
					output_line += '\n'
					vcf_output += output_line
				o = file(outfile,"w")
				o.write(vcf_output)
				o.close()
				print "\n... done"
				if add_vcf:
					return outfile

		except:
			print "\nCouldn't print VCF, check option -v has been set correctly"

	def mergeVCF(vcf_to_add, pre):
		print "\nCreating merged VCF with PASS/FAIL information..."
		snp_filter = "PASS"
		vcf_out = []
		vcf_out = getVCF(vcf_out, vcf_to_add[-1], snp_filter)
		for i in range(len(vcf_to_add)-1):
			snp_filter = "N:R"+str(len(vcf_to_add)-(i+1))
			vcf_out = getVCF(vcf_out, vcf_to_add[len(vcf_to_add)-(i+2)], snp_filter)
		o = file(pre + '_gnr.vcf',"w")
		for line in vcf_out:
			o.write(line)
		o.close()
		print "\n... done"

	def getVCF(vcf_in, vcf_new, snp_filter):
		new_vcf = file(vcf_new, "r")
		if vcf_in == []:
			vcf_out = new_vcf.readlines()
			new_vcf.close()
		else:
			vcf_add = new_vcf.readlines()
			new_vcf.close()
			vcf_out = []
			header_count = 0
			for line in vcf_in:
				if line.startswith("#"):
					header_count += 1
			for i in range(header_count):
				if i != 1 and i != header_count-1: 
					vcf_out.append(vcf_in[i])
				elif i == 1:##source
					new_vcf_out = vcf_in[i][:-1]
					source = vcf_add[i].split('=')
					new_vcf_out += ' + '+source[-1]
					vcf_out.append(new_vcf_out)
				else:
					new_vcf_out = '##FILTER=<ID='+snp_filter.lstrip('N:')+',Description="SNP filtered out in round '+snp_filter.lstrip('N:R')+'">\n'
					vcf_out.append(new_vcf_out)
					vcf_out.append(vcf_in[i])
			i = header_count
			j = header_count
			while i < len(vcf_in) and j < len(vcf_add):
				snp_in = vcf_in[i].split('\t')
				snp_add = vcf_add[j].split('\t')
				if int(snp_add[1]) == int(snp_in[1]):
					vcf_out.append(vcf_in[i])
					i += 1
					j += 1
				elif int(snp_add[1]) > int(snp_in[1]): #this shouldn't occur, but just in case
					vcf_out.append(vcf_in[i])
					i += 1
				elif int(snp_add[1]) < int(snp_in[1]):
					for k in range(len(snp_add)):
						if k != 6 and k != 0:
							new_vcf_out += '\t' + snp_add[k]
						elif k != 0:
							new_vcf_out += '\t' + snp_filter
						else: 
							new_vcf_out = snp_add[k]
					vcf_out.append(new_vcf_out)
					j += 1
			if i >= len(vcf_in):
				while j < len(vcf_add):
					snp_add = vcf_add[j].split('\t')
					for k in range(len(snp_add)):
						if k != 6 and k != 0:
							new_vcf_out += '\t' + snp_add[k]
						elif k != 0:
							new_vcf_out += '\t' + snp_filter
						else: 
							new_vcf_out = snp_add[k]
					vcf_out.append(new_vcf_out)
					j += 1
			elif j >= len(vcf_add): # Again, shouldn't happen, but just in case
				while i < len(vcf_in):
					vcf_out.append(vcf_in[i])
					i += 1
		return vcf_out

	def filterCore(snptable, strains, pre, options, outgroups):
		if options.refseq=="":
			print "\nNo reference genbank file specified (-r), can't do core SNP filtering"
		else:
			try:
				print "\nFiltering SNPs based on genes in core genome in core strains..."
				core_strains_file_name = options.core_strains
				core_strains = []
				if core_strains_file_name == "":
					for strain in strains:
						if strain not in outgroups:
							core_strains.append(strain)
				else:
					core_strains_file = open(core_strains_file_name, 'r')
					for line in core_strains_file:
						core_strain = line.rstrip('\n')
						if core_strain not in outgroups and core_strain in strains:
							if core_strain not in core_strains:
								core_strains.append(core_strain)
					core_strains_file.close()
				print "\nReading gene features from reference " + options.refseq
				## READ IN GENBANK FILE 
				Ref_Passed = True
				handle = open(options.refseq,"r")
				if options.queryseq=="":
					try:
						record = SeqIO.read(handle, "genbank")
						sequence = record.seq
						geneannot = record.features
						mapped = record.name
					except:
						Ref_Passed = False
						print "\nCheck reference sequence for multiple records: can't do core SNP filtering"
				else:
					records = SeqIO.parse(handle, "genbank")
					Ref_Passed = False
					mapped = options.queryseq
					for item in records:
						if item.name == mapped:
							record = item
							sequence = SeqRecord(item.seq)
							geneannot = item.features
							Ref_Passed = True
					if Ref_Passed == False:
						print "\nCheck reference sequence: queryseq (-q) not found"
				handle.close()
				if Ref_Passed==True:
					gene_list = []
					gene_position = []
					for f in geneannot:
						if f.type == "CDS":
							start = f.location.nofuzzy_start
							stop = f.location.nofuzzy_end
							sysid = mapped+";"+str(start+1)+'-'+str(stop+1)
							f.qualifiers['sysid'] = [sysid]
							if 'locus_tag' in f.qualifiers:
								locus_tag = f.qualifiers['locus_tag'][0]
							else:
								#if the locus_tag is missing from the genbank record make up a tag as RedDog does
								locus_tag = "tag_" + str(start)+'-'+str(stop)
							gene_list.append(locus_tag)
							gene_position.append([start, stop])
					if options.gene_coverage == "":
						print "\nNo gene coverage file specified (-z), can't do core SNP filtering"
					else:						
						core_coverage = float(options.core_coverage)
						if core_coverage < 0 or core_coverage > 1:
							print "\nCore coverage (-Z) outside expect range (0 to 1), can't do core SNP filtering"
						else:
							core_strain_index = []
							ordered_core_genes = []
							gene_coverage = open(options.gene_coverage, 'r')
							for line in gene_coverage:
								# header
								line = line.rstrip('\n')
								if line.startswith("replicon__gene"):
									mapped_strains = line.split(',')
									for i in range(1,len(mapped_strains)):
										if mapped_strains[i] in core_strains:
											core_strain_index.append(i)
									core_count_test = len(core_strain_index)
								#core genes
								if line.startswith(mapped):
									entry = line.split(',')
									core_count = 0
									for i in core_strain_index:
										if float(entry[i])/100 >= core_coverage:
											core_count += 1
									if core_count == core_count_test:
										tag = entry[0].split('__')
										locus_tag = tag[1]
										i = gene_list.index(locus_tag)
										start = min(gene_position[i][0],gene_position[i][1])
										stop = max(gene_position[i][0],gene_position[i][1])
										if start == gene_position[i][1]: #ie if gene in reverse position
											start += 1
											stop += 1
										ordered_core_genes.append([start,stop])
							gene_coverage.close()
							ordered_core_genes.sort(key=operator.itemgetter(1))
							ordered_core_genes.sort(key=operator.itemgetter(0))
							ordered_snp_list = []
							for snp in range(len(snptable)): # cycle over SNPs
								ordered_snp_list.append([int(snptable[snp][0]),snp])
							ordered_snp_list.sort(key=operator.itemgetter(0))

							to_remove = []
							i = 0
							j = 0
							while i < len(ordered_snp_list) and j < len(ordered_core_genes):
								if ordered_snp_list[i][0] < ordered_core_genes[j][0]:
									to_remove.append(ordered_snp_list[i][1])
									i += 1
								elif ordered_snp_list[i][0] < ordered_core_genes[j][1]:
									i += 1
								else:
									j += 1
							if i < len(ordered_snp_list):
								while i < len(ordered_snp_list):
									to_remove.append(ordered_snp_list[i][1])
									i += 1
							if to_remove != []:
								filtered_snptable = []
								for snp in range(len(snptable)):
									if snp not in to_remove:
										filtered_snptable.append(snptable[snp])
								snptable = filtered_snptable
							pre += "_core" + str(core_coverage)
							outfile = pre + ".csv"
							o = open(outfile,"w")	
							o.write(",".join(["Pos"] + strains) + "\n") # write header
							for snp in range(len(snptable)):
								o.write(snptable[snp][0])
								for strain in range(len(strains)):
									o.write(","+snptable[snp][1][strain])
								o.write("\n")
							o.close()
							print "\n... " + str(len(to_remove)) + " SNPs removed as not in core genes"
							print "... " + str(len(snptable)) + " SNPs passed filter; printed to " + outfile
			except:
				print "\nCouldn't filter SNPs based on genes in core genome, check the following option(s):"
				if options.core_strains != "":
					print "    core_strains  -L " + options.core_strains
				print "    gene_coverage  -z " + options.gene_coverage
				print "    core_coverage -Z " + options.core_coverage
		return pre, snptable

	# run module and return updated values for snptable, strains and pre
	# aln,fasttree,filter,rax,coding
	def runModule(m, snptable, strains, pre, options, complement, outgroups, add_vcf, vcf_to_add):
		if m == "aln":
			# print mfasta alignment of the current table
			printFasta(snptable, strains, pre + ".mfasta")
		elif m == "cons":
			pre, snptable = filterCons(snptable, strains, pre, options, outgroups)
		elif m == "core":
			pre, snptable = filterCore(snptable, strains, pre, options, outgroups)
		elif m == "fasttree":
			runFasttree(pre, snptable, strains)
		elif m == "filter":
			pre, snptable = filter(snptable, strains, pre, options) # return filtered snp table
		elif m == "clean":
			pre, snptable = cleanSNPs(snptable, strains, pre, options)
		elif m == "vcf":
			if add_vcf:
				vcf = printVCF(snptable, strains, pre, options, add_vcf)
				if vcf != "":	
					vcf_to_add.append(vcf)
				else:
					print "\nVCF not added: check -v option"
			else:
				printVCF(snptable, strains, pre, options, add_vcf)
		elif m == "rax":
			runRax(pre, options, snptable, strains)
		elif m == "coding":
			runCoding(pre, snptable, options, complement)
		return pre, snptable, vcf_to_add
		
	### MAIN PROCESS
		
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-': '-'} 
	
	# set up variables
	(dir,filename) = os.path.split(options.snptable)
	(pre,ext) = os.path.splitext(filename) # pre = current file prefix, updated if file is filtered to exclude SNPs
	
	pre = options.prefix + pre # add user specified prefix

	if options.directory != "":
		if options.directory[:-1] != '/':		
			pre = options.directory + '/' + pre
		else:
			pre = options.directory + pre

	outgroups = [] # list of outgroups provided
	if options.outgroup != "":
		outgroups = options.outgroup.split(",")

	# read in raw SNP table
	snptable, strains, pre = readSNPTable(options.snptable,options.outgroup,options.subset,pre)

	#create vcf of raw SNP table if needed
	vcf_to_add = []
	add_vcf = options.add_vcf
	if add_vcf:
		vcf = printVCF(snptable, strains, pre, options, add_vcf)
		if vcf != "":	
			vcf_to_add.append(vcf)
		else:
			print "\nFirst VCF not added; 'add_vcf' set to 'False': check -v option"
			add_vcf = False
	
	# run modules in order
	for m in options.modules.split(","):
		pre, snptable, vcf_to_add = runModule(m, snptable, strains, pre, options, complement, outgroups, add_vcf, vcf_to_add)

	#create combined VCF of with PASS/FAIL if needed
	if add_vcf and len(vcf_to_add) == 1:
		print "\nOnly one VCF added: ensure at least one 'vcf' module in -m after a filtering stage (filter, cons or clean)"
	if add_vcf and len(vcf_to_add) == 0:
		print "\nNo VCF added: check option reference name -v"
	elif add_vcf:
		mergeVCF(vcf_to_add, pre)		
	print ""
#	print ('Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
