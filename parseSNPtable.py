# read in SNP alleles (csv format)
# and a list of coordinates of regions to exclude (genbank, gff or 2-column CSV table format)
# output list of SNP loci and their alleles
# optionally launch fasttree or raxml

import os, sys, subprocess, string, re, random
import collections
from optparse import OptionParser
from datat import Datat
from datat import load_csv
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

	# modules to run
	parser.add_option("-m", "--modules", action="store", dest="modules", help="modules to run, comma separated list (aln,filter,coding,phy,rax,fasttree,raxRaw,fasttreeRaw)", default="aln,filter,coding,phy,rax,fasttree")

	# snp filtering
	parser.add_option("-s", "--snptable", action="store", dest="snptable", help="SNP table (CSV)", default="")
	parser.add_option("-x", "--exclude", action="store", dest="excluded", help="file of excluded regions", default="")
	
	# coding consequences
	parser.add_option("-r", "--refseq", action="store", dest="refseq", help="reference sequence file (gbk)", default="")
	parser.add_option("-f", "--genefeatures", action="store", dest="genefeatures", help="feature types for protein coding genes (default CDS; can be multiple comma-sep)", default="CDS")
	parser.add_option("-e", "--excludefeatures", action="store", dest="excludefeatures", help="feature types to exclude (default gene,misc_feature)", default="gene,misc_feature")
	parser.add_option("-i", "--identifier", action="store", dest="identifier", help="unique identifier for features (locus_tag)", default="locus_tag")

	# phylogeny
	parser.add_option("-F", "--fasttree", action="store", dest="fasttree", help="run fasttree (default 0)", default="0")
	parser.add_option("-R", "--raxml", action="store", dest="raxml", help="run RAxML (default 0)", default="0")

	return parser.parse_args()


if __name__ == "__main__":

	(options, args) = main()
	
	# helper functions
	
	def table2mfasta(table):
		(pre,ext) = os.path.splitext(table)
		f=file(table,"r")
		l=f.readline()
		snpcol = l.rstrip().split(",")[0]
		snptable = load_csv(table,namescol=snpcol) # snptable object
		m = open(pre + ".mfasta","w")
		for strain in snptable.columns:
			if strain != "Allele":
				m.write(">" + strain + "\n")
				seq = ""
				for snp in snptable: # cycle over SNPs
					seq = seq + str(snptable[snp][strain])
				m.write(seq + "\n")
		m.close()
		print "\nAlignment of SNPs in table " + table + " written to " + pre + ".mfasta"
		return(snptable)
		
	def getCodons(genestart,genestop,genestrand,snp,derived,ancestral):
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
		
	def runFasttree(mfasta):
		(pre,ext) = os.path.splitext(mfasta)
		cmd = "#!/bin/bash\n#PBS -N ft" + pre
		cmd += "\n#PBS -l walltime=12:00:00"
		cmd += "\n#PBS -l pvmem=4gb"
		cmd += "\ncd " + os.getcwd()
		cmd += "\nmodule load fasttree\n"
		cmd += "FastTree -gtr -gamma -nt " + mfasta + " > " + pre + ".tree"
		print "\nRunning fasttree on " + mfasta + ", using job script:"
		print cmd
		os.system('echo "' + cmd + '" | qsub')
		
	def runRax(mfasta,cutoff,numreps,N):
		print "\nFiltering alignment " + mfasta + " to remove SNPs with <" + str(cutoff*100) + "% coverage and converting to phylip format for RAxML"
		(pre,ext) = os.path.splitext(mfasta)
		filteredFilePhylip = pre + ".coverageFiltered60.phy"
		input_handle = open(mfasta, "rU")
		output_handle = open(filteredFilePhylip, "w")
		alignments = AlignIO.parse(input_handle,"fasta",alphabet=generic_dna)
		filtered = []
		for a in alignments:
			for i in a:
				sequence = list(i.seq)
				bases = len(sequence)
				missing = sequence.count('-')
				covered = float(bases-missing)/bases
				if covered >= cutoff:
					filtered.append(i)
		filteredalignments = MultipleSeqAlignment(filtered)
		AlignIO.write(filteredalignments, output_handle, "phylip")
		output_handle.close()
		input_handle.close()
	
		print "\nRunning RAxML on " + filteredFilePhylip + ", using " + str(numreps) + " replicate runs with " + str(N) + " bootstraps each"
		for rep in range(1,numreps):
			seed = random.randint(10000,1000000)
			if seed % 2 == 0:
				seed += 1
	
			cmd = "#!/bin/bash\n#PBS -N rax" + pre + str(rep)
			cmd += "\n#PBS -l walltime=48:00:00"
			cmd += "\n#PBS -l pvmem=4gb"
			cmd += "\n#PBS -l procs=8"
			cmd += "\ncd " + os.getcwd()
	
			# load modules
			cmd += "\nmodule load raxml-gcc/7.2.8.alpha\n"
	
			# commands
			cmd += "mpiexec raxmlHPC-MPI -s " + filteredFilePhylip
			cmd += " -n " + pre + "." + str(rep) + " -f a -m GTRGAMMA -x " + str(seed) + " -N " + str(N)
	
			os.system('echo "' + cmd + '" | qsub')



	# modules to run
	modules = options.modules.split(",")
	
	# set up variables
	snptable = None # current SNP table
	(dir,filename) = os.path.split(options.snptable)
	(pre,ext) = os.path.splitext(filename) # pre = current file prefix, updated if file is filtered to exclude SNPs
	
	
	# raw SNP set
	
	if "aln" in modules:
		# generate alignment of unfiltered SNPs
		snptable = table2mfasta(options.snptable)
		
	if "fasttreeRaw" in modules or ("fasttree" in modules and "filter" not in modules):
		if not os.path.exists(pre + ".mfasta"):
			snptable = table2mfasta(options.snptable) # generate mfasta
		runFasttree(pre + ".mfasta") # run fasttree
	
	if "raxRaw" in modules or ("rax" in modules and "filter" not in modules):
		if not os.path.exists(pre + ".mfasta"):
			snptable = table2mfasta(options.snptable) # generate mfasta
		runRax(pre + ".mfasta",0.6,10,1000) # run RAxML
			
	

	if "filter" in modules:
	
		# filter SNPs from excluded regions
		pre = pre + "_filteredAlleles"
		
		# parse excluded genomic regions
		if options.excluded.split('.')[-1] == "gbk" or options.excluded.split('.')[-1] == "gb":
			excluded = []
			handle = open(options.excluded,"r")
			record = SeqIO.read(handle, "genbank")
			for region in record.features:
				for i in range(int(region.location.nofuzzy_start),int(region.location.nofuzzy_end)+1):
					excluded.append(i)
		else: 
			start = 0
			stop = 1
			delim = "," # assume csv
			if options.excluded.split('.')[-1] == "gff":
				start = 3
				stop = 4
				delim = "\t" # gff is tab-delim
			e = file(options.excluded, "r")
			excluded = {}
			for line in e:
				fields = line.rstrip().split(delim)
				for i in range(int(fields[start]),int(fields[stop])+1):
					excluded[i] = 1
			e.close()
	
		# read in snp table & print lines that are not excluded
		bases = ['A','C','G','T']
		print "\nFiltering SNPs"
		snpcount = 0
		o = open(pre + ".csv","w")
		s = file(options.snptable, "r")
		header = 1
		for line in s:
			if header > 0:
				header = 0
				o.write(line) # write header row = strain names
			else:
				fields = line.split(',')
				pos = fields.pop(0)
				base_count = 0
				for i in bases:
					if fields.count(i) > 0:
						base_count += 1
				if base_count >= 2 and int(pos) not in excluded:
					o.write(line)
					snpcount += 1
		s.close()
		o.close()
		print "  " + str(snpcount) + " SNPs passed filter; printed to " + pre + ".csv"
		
		# read in filtered table and print as mfasta alignment
		snptable = table2mfasta(pre + ".csv")
		
		
		# launch phylogenetic analysis
		
		if "fasttree" in modules:
			runFasttree(pre + ".mfasta")
			
		if "rax" in modules:
			runRax(pre + ".mfasta",0.6,10,1000)
		
	### end filtering module
	

	
	if "coding" in modules:
	
		if options.refseq=="":
			DoError("No reference genbank file specified")	
			
		complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-': '-'} 
		
		genefeatures = options.genefeatures.split(",")
		excludefeatures = options.excludefeatures.split(",")
			
		if snptable is None:
			# read in the input SNP file; otherwise if alignment or filtering has been requested the appropriate file is already loaded
			snptable = load_csv(options.snptable,namescol="Pos")
			
		# order SNPs
		snp_list_ordered = []
		for snp in snptable:
			snp_list_ordered.append(int(snp))
		snp_list_ordered.sort()
		
		print "\nReading gene features from reference " + options.refseq
		
		# check coding consequences and generate genbank file of SNP loci
	
		## READ IN GENBANK FILE 
		handle = open(options.refseq,"r")
		record = SeqIO.read(handle, "genbank")
		sequence = record.seq
		geneannot = record.features
			
		print "\nDetermining coding changes"
			
		## GET CONSEQUENCES FOR SNPS and WRITE SNP ANNOTATION FILE
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
			if "-" in allele_list:
				allele_list.remove("-")
			if ref_allele in allele_list:
				allele_list.remove(ref_allele)
			if len(allele_list)>0:
				for alt_allele in allele_list:
					hit = 0 # initialize
					for feature in geneannot:
						if int(snp) in feature and feature.type != "source" and feature.type not in excludefeatures:
							hit = 1
							start = int(feature.location.nofuzzy_start) # feature start
							stop = int(feature.location.nofuzzy_end) # feature stop
							id = ""
							product = ""
							if options.identifier in feature.qualifiers:
								id = feature.qualifiers[options.identifier][0]
							if 'product' in feature.qualifiers:
								product = feature.qualifiers['product'][0]
							if feature.type in genefeatures:
								# get coding effect of coding features
								(ancestral_codon,derived_codon,ancestralAA,derivedAA,posingene,posincodon)=getCodons(start,stop,feature.strand,int(snp),alt_allele,ref_allele)
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
								if feature.strand == 1:
									posingene = snp-start+1
								else:
									posingene = stop-snp+1
								o.write("\t".join([str(snp),str(ref_allele),str(alt_allele),feature.type,id,"","","","",product,str(posingene),"","","\n"]))
								record.features.append(SeqFeature(FeatureLocation(int(snp)-1,int(snp)), type="variation", strand=1, qualifiers = {'note' : ["SNP " + ref_allele + "->" + alt_allele + " in non-CDS feature" ]}))
					if hit == 0:
						# SNP is intergenic
						intergenic_count += 1
						o.write("\t".join([str(snp),str(ref_allele),str(alt_allele),"intergenic","","","","","","","","\n"]))
						record.features.append(SeqFeature(FeatureLocation(int(snp)-1,int(snp)), type="variation", strand=1, qualifiers = {'note' : ["intergenic SNP " + ref_allele + "->" + alt_allele]}))
		o.close()
		SeqIO.write(record, pre + ".gbk", "genbank")
		
		print "  " + str(ns_count) + " nonsyonymous, " + str(syn_count) + " synonymous, " + str(other_feature_count) + " in other features, " + str(intergenic_count) + " in non-coding regions"
		print "  Coding consequences written to file " + pre + "_consequences.txt"
		print "  SNP loci annotated in genbank file " + pre + ".gbk"
