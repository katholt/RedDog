#!/usr/bin/env python
#
# Read in SNP alleles (csv format)
#  can take a file containing a list of strains to include (-l), otherwise all are included
#  can take a list of outgroups (-o) => if specified SNPs that are non-variable in in the ingroup are removed on read-in; and variation in the outgroups is ignored in assessing conservation
#  Note it doesn't matter whether or not the outgroup is included in the list of strains.
# The table is then parsed as specified by the modules in -m, these include:
#   aln - convert to fasta alignment
#   filter - filter SNPs that are included/excluded in regions specified via -x (genbank, gff or 2-column CSV table format)
#   cons - filter SNP positions that are not conserved above a cutoff specified via -c (e.g. -c 0.99 -> all snps with >1% missing alleles is filtered out)
#   fasttree - submit a fasttree job to SLURM
#   rax - submit threaded RAxML jobs to SLURM (can specify walltime, memory, threads, number of boostraps and number of replicate runs)
# Any number of these modules can be supplied in any order; the order they are given is the order they will be run
#   specify modules in a comma-separated list, e.g. '-m filter,cons,aln,rax' will run region filtering, then conservation filter, then make a fasta alignment and use this as input to 5 raxml jobs with 100 bootstraps using 8 threads
#
# This version also handles multiple sequence genbank files for coding entries. The sequence in the genbank relevant to the SNP table must be specified (-q queryseq).  
# It is also quicker generating coding consequences. 

# NOTE this script does not submit fasttree or RAxML jobs to SLURM
#
# Authors - Kat holt (kholt@unimelb.edu.au)
#         - David Edwards (d.edwards2@student.unimelb.edu.au) 
#
# Example command on barcoo:
'''
module load python-gcc/2.7.5
python /vlsci/VR0082/shared/code/holtlab/parseSNPtable_multigbk.py -s snps.csv -p prefix -r genbank -q queryseq -m aln, coding, rax
'''
#
# Last modified - Mar 25, 2014
# Changes:
#	 15/10/13 - added strain subset option
#    25/03/14 - added multiple sequence genbank file handling and improved 'coding' option performance

import os, sys, subprocess, string, re, random
import collections
from optparse import OptionParser
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)


	parser.add_option("-p", "--prefix", action="store", dest="prefix", help="prefix to add to output files (default none)", default="")
	parser.add_option("-d", "--directory", action="store", dest="directory", help="directory to send output files (default none)", default="")
	
	# modules to run
	parser.add_option("-m", "--modules", action="store", dest="modules", help="modules to run, comma separated list in order (filter,cons,aln,coding)", default="filter,cons,aln,coding")

	# snp filtering
	parser.add_option("-s", "--snptable", action="store", dest="snptable", help="SNP table (CSV)", default="")
	parser.add_option("-x", "--regions", action="store", dest="regions", help="file of regions to include/exclude (gbk)", default="")
	parser.add_option("-y", "--include", action="store", dest="include", help="include (default exclude)", default="exclude")
	parser.add_option("-g", "--gapchar", action="store", dest="gapchar", help="gap character (default -, could be N)", default="-")
	parser.add_option("-c", "--conservation", action="store", dest="conservation", help="minimum conservation across samples required to retain SNP locus (default 0.99)", default="0.99")
	parser.add_option("-o", "--outgroup", action="store", dest="outgroup", help="comma separated list; outgroup strains (alleles will be included but not sites that vary only in outgroups)", default="")
	parser.add_option("-l", "--subset", action="store", dest="subset", help="file containing list of strains to include (one per line), otherwise all strains included", default="")
	
	# coding consequences
	parser.add_option("-r", "--refseq", action="store", dest="refseq", help="reference sequence file (gbk)", default="")
	parser.add_option("-q", "--queryseq", action="store", dest="queryseq", help="query sequence in reference sequence file (multisequence gbk)", default="")
	parser.add_option("-f", "--genefeatures", action="store", dest="genefeatures", help="feature types for protein coding genes (default CDS; can be multiple comma-sep)", default="CDS")
	parser.add_option("-e", "--excludefeatures", action="store", dest="excludefeatures", help="feature types to exclude (default gene,misc_feature)", default="gene,misc_feature")
	parser.add_option("-i", "--identifier", action="store", dest="identifier", help="unique identifier for features (locus_tag)", default="locus_tag")


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
		
		snptable = collections.defaultdict(dict) # key 1 = snp, key 2 = strain
		strains = [] # strains from header
		ignored = []
		pre += "_var"
		f = file(infile, "r")
		o = file(pre + ".csv","w")
		o.write("Pos")
		for line in f:
			fields = line.rstrip().split(',')
			if len(strains)==0:
				strains = fields
				if len(strainlist) == 0:
					for i in range(1,len(strains)):
						strain = strains[i]
						if strain not in outgroups:
							strainlist.append(strain) # retain all strains (except outgroups)
				# remove strains from the strainlist if we have not encountered them in the actual table
				for strain in strainlist:
					if strain not in strains:
						strainlist.remove(strain)
				# print header for new table
				for i in range(1, len(fields)):
					strain = strains[i]
					if strain in strainlist or strain in outgroups:
						o.write(","+strain)
						
				o.write("\n")
			else:
				snp = fields[0]
				alleles_ingroup = [] # alleles for this snp
				alleles_outgroup = [] # alleles for outgroups
				for i in range(1,len(fields)):
					strain = strains[i] # col name
					if strain in strainlist and strain not in outgroups:
						alleles_ingroup.append(fields[i]) # add this allele to the list
					elif strain in outgroups:
						alleles_outgroup.append(fields[i]) # add this allele to the outgroup list
						if strain not in outgroups_used:
							outgroups_used.append(strain)
				if isVar(alleles_ingroup):
					o.write(str(snp))
					# variable in ingroup, store alleles for included strains
					for i in range(1, len(fields)):
						strain = strains[i]
						if strain in strainlist or strain in outgroups:
							snptable[snp][strain] = fields[i].upper()
							o.write(","+fields[i].upper())
					o.write("\n")
				else:
					ignored.append(snp)
		f.close()
					
		strains.pop(0) # remove SNP column header

		print "\n... finished reading " + str(len(snptable)) + " variable SNPs in " + str(len(strainlist)) + " ingroup strains"
		print "... ignoring " + str(len(ignored)) + " SNPs that are non-variable among these ingroup strains"
		
		return(snptable, strainlist + outgroups_used, pre) # include outgroups that appear in the snptable in strainlist
	
	def printFasta(snptable, strains, outfile):
		print "\nPrinting alignment to file " + outfile
		o = file(outfile,"w")
		for strain in strains:
			o.write(">" + strain + "\n")
			seq = ""
			for snp in snptable: # cycle over SNPs
				seq = seq + str(snptable[snp][strain])
			o.write( seq + "\n")
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
			for snp in snptable:
				snp_list_ordered.append(int(snp))
			snp_list_ordered.sort()		
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
					for slice in range((len(sequence)/slice_size)+1):
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
					ref_allele = sequence[int(snp)-1]
					allele_list = []
					for strain in snptable[str(snp)]:
						if snptable[str(snp)][strain] not in allele_list:
							allele_list.append(snptable[str(snp)][strain])
					if options.gapchar in allele_list:
						allele_list.remove(options.gapchar)
					if ref_allele in allele_list:
						allele_list.remove(ref_allele)
					if len(allele_list)>0:
						for alt_allele in allele_list:
							hit = 0 # initialize
							snp_slice = int(snp)/slice_size
							if feature_slice[snp_slice] != []:
								for feature_index in feature_slice[snp_slice]:
									if int(snp) in geneannot[feature_index[2]]:
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
											(ancestral_codon,derived_codon,ancestralAA,derivedAA,posingene,posincodon)=getCodons(start,stop,geneannot[feature_index[2]].strand,int(snp),alt_allele,ref_allele,sequence,complement)
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
											record.features.append(SeqFeature(FeatureLocation(int(snp)-1,int(snp)), type="variation", strand=1, qualifiers = {'note' : [note]}))
											o.write("\t".join([str(snp),ref_allele,alt_allele,change,id,str(ancestral_codon),str(derived_codon),str(ancestralAA),str(derivedAA),product,str(posingene),str(codon_number),str(posincodon),"\n"]))
										else:
											# non-protein coding feature
											other_feature_count += 1
											if geneannot[feature_index[2]].strand == 1:
												posingene = snp-start+1
											else:
												posingene = stop-snp+1
											o.write("\t".join([str(snp),str(ref_allele),str(alt_allele),geneannot[feature_index[2]].type,id,"","","","",product,str(posingene),"","","\n"]))
											record.features.append(SeqFeature(FeatureLocation(int(snp)-1,int(snp)), type="variation", strand=1, qualifiers = {'note' : ["SNP " + ref_allele + "->" + alt_allele + " in non-CDS feature" ]}))
							if hit == 0:
								# SNP is intergenic
								intergenic_count += 1
								o.write("\t".join([str(snp),str(ref_allele),str(alt_allele),"intergenic","","","","","","","","\n"]))
								record.features.append(SeqFeature(FeatureLocation(int(snp)-1,int(snp)), type="variation", strand=1, qualifiers = {'note' : ["intergenic SNP " + ref_allele + "->" + alt_allele]}))
				o.close()
				SeqIO.write(record, pre + ".gbk", "genbank")
				print "\n... " + str(ns_count) + " nonsyonymous, " + str(syn_count) + " synonymous, " + str(other_feature_count) + " in other features, " + str(intergenic_count) + " in non-coding regions"
				print "... coding consequences written to file " + pre + "_consequences.txt"
				print "... SNP loci annotated in genbank file " + pre + ".gbk"

			
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
					for i in range(int(region.location.nofuzzy_start),int(region.location.nofuzzy_end)+1):
						regions.append(i)
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
					for i in range(int(fields[start]),int(fields[stop])+1):
						if i not in regions:
							regions.append(i)
				e.close()
			total_masked_bases = len(regions)
			# filter snps
			o = file(pre + "_regionFiltered.csv","w")
			o.write(",".join(["Pos"] + strainlist) + "\n") # write header
			print "\nFiltering SNPs that are ",
			if options.include != "include":
				# excluding snps 
				print "located in excluded regions",
			else:
				print "located outside included regions",
			print "totalling " + str(total_masked_bases) + " bases",
			print "specified in file " + options.regions
			snpcount = 0
			to_remove = [] # list of snps to remove
			for snp in snptable:
				keep = False
				if options.include != "include" and int(snp) not in regions:
					keep = True
				elif options.include == "include" and int(snp) in regions:
					keep = True
				if keep:
					alleles = snptable[snp]
					o.write(snp)
					for strain in strainlist:
						o.write(","+alleles[strain])
					o.write("\n")
					snpcount += 1
				else:
					to_remove.append(snp)
			o.close()
			pre += "_regionFiltered"
			print "... " + str(snpcount) + " SNPs passed filter; printed to " + pre + ".csv"
			for snp in to_remove:
				del snptable[snp]
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
			n = len(strainlist) # total strains
			to_remove = [] # list of snps to remove
			snp_list_ordered = []
			for snp in snptable:
				snp_list_ordered.append(int(snp))
			snp_list_ordered.sort()		
			snpcount = 0
			for snp in snp_list_ordered:
				alleles = snptable[str(snp)] # dictionary
				if len(outgroups) == 0:
					allele_list = alleles.values() # all alleles
				else:
					allele_list = []
					for strain in alleles:
						if strain not in outgroups:
							allele_list.append(alleles[strain])
				numgaps = allele_list.count(options.gapchar) # number of gaps at this position, excluding outgroups
				callrate = 1 - float(numgaps) / n
				if callrate >= float(options.conservation):
					o.write(str(snp))
					for strain in strainlist:
						o.write(","+alleles[strain])
					o.write("\n")
					snpcount += 1
				else:
					to_remove.append(snp) # remove SNP from table
			o.close()
			print "\n... " + str(snpcount) + " SNPs passed filter; printed to " + pre + ".csv"
			for snp in to_remove:
				del snptable[str(snp)]
		except:
			print "\nCouldn't filter SNPs based on missing alleles, couldn't understand the proportion given: -c " + options.conservation
		return pre, snptable
		
	# run module and return updated values for snptable, strains and pre
	# aln,fasttree,filter,rax,coding
	def runModule(m, snptable, strains, pre, options, complement, outgroups):
		if m == "aln":
			# print mfasta alignment of the current table
			printFasta(snptable, strains, pre + ".mfasta")
		elif m == "cons":
			pre, snptable = filterCons(snptable, strains, pre, options, outgroups)
		elif m == "filter":
			pre, snptable = filter(snptable, strains, pre, options) # return filtered snp table
		elif m == "coding":
			runCoding(pre, snptable, options, complement)
		return pre, snptable
		
	### MAIN PROCESS
		
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-': '-'} 
	
	# set up variables; pre = current prefix; snptable = current snptable dictionary of dictionaries
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
	
	# run modules in order
	for m in options.modules.split(","):
		pre, snptable = runModule(m, snptable, strains, pre, options, complement, outgroups)
		
	print ""
