'''
Configuration file for RedDog.py V1beta.10.4
-------------------------------

Copyright (c) 2016 David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)

Essential pipeline variables.
'''
no_check = True
reference = ""
sequences = ""
output = ""
out_merge_target = ""

'''
force_tree and force_no_tree
The pipeline can produce a FastTree. If there are more than 500 isolates, the tree generation
will be switched off. If you want a tree for larger data sets, set force_tree to 'True'.
If you want to turnoff the tree entirely, set force_no_tree to 'True'.
(this overrides force_tree)
'''
force_tree = False
#force_tree = True

force_no_tree = False
#force_no_tree = True

'''
Notes:

'no_check' is for switching off the user check during the start of the run. This enables 
the pipeline to be run as a single job (with lots of cpus) on a distributed system

'reference' and 'sequences'

Reference can be GenBank or fasta format - if GenBank format, this will be converted
to a fasta version for mapping.

If the GenBank file is not given, the gene cover and depth matrices for genes
will not be generated and nor will SNP consequences.

If you don't have the GenBank record, or don't want the above matrices
to be generated, enter a fasta format reference instead.

'''
#Test Sets
#reference = "/full_path_to/pipeline_test_sets/reference/NC_007384.gbk"
#reference = "/full_path_to/pipeline_test_sets/reference/NC_007384_with_plasmid.gbk"
#reference = "/full_path_to/pipeline_test_sets/reference/NC_007384_with_plasmid.fasta"
#sequences = "/full_path_to/pipeline_test_sets/*.fastq.gz"
#sequences = "/full_path_to/pipeline_test_sets/extra/*.fastq.gz"

# You can now also combine sequences from different folders into the same run...
#sequences = ["/full_path_to/pipeline_test_sets/*.fastq.gz", "/full_path_to/pipeline_test_sets/extra/*.fastq.gz"]

'''
'output' directory:
full path name including final "/"

For large data sets run the output to the scratch disk area and save the final output to
your shared directory (or contagion if you have access)
e.g. output = "/scratch/VR0082/a_folder/<ref>_<version>_<date>/"
'''
#output = "/full_path_to/<your_directory>/RedDog_output/<ref>_<version>_<date>/"

'''
Directory to merge output with ('out_merge_target'):

When running new analysis set to null string.

Otherwise set to the directory you want to merge with.

You can only merge a prior run with a new run (not two prior runs)
This merge target folder must have the bams and indexes in one
sub-folder (/bam) and the vcfs in another (/vcf). There also must
be a sequence_list.txt file.

The 'output' folder (see above) for a merge run should NOT exist prior to the run,
and will be deleted at completion of the pipeline.

Set to empty string for no merging (i.e. new run).

'''
#out_merge_target = ""
#out_merge_target = "/full_path_to/<your_directory>/RedDog_output/<ref>_<version>_<date>/"

'''
If none of the following are set or changed, the default settings will be used.
###############################################################################

Choose type of input sequences:    IT for ion torrent (single reads),
                                   PE for Illumina pair-end reads
                                or SE for Illumina single-end reads

readType = "PE" or "SE" or "IT"

Each read type has a particular pattern you need to follow for use in the pipeline:

SE: *.fastq.gz                      (i.e. must be gzipped)
PE: *_1.fastq.gz and *_2.fastq.gz   (i.e. must be gzipped with forward and reverse in separate files)
IT: *_in.iontor.fastq.gz            (i.e. must be gzipped)

'''
readType = "PE"
#readType = "SE"
#readType = "IT"

'''
Run Type - if set to a null string, the number of replicons in the reference will determine
the run type:
    1 - 100 replicons (e.g. reference genome + plasmids + phage)   - phylogeny run type
    > 100 replicons (e.g. multifasta pangenome)                    - pangenome run type

The user can override the run type, by setting it below. The run types are described in
more detail in the instructions.

'''
runType = ""
#runType = "pangenome"
#runType = "phylogeny"

'''
For a pangenome run, the SNPs will only be called for the largest replicon - this assumes
the core genome is in this replicon. The user can define an alternative replicon
(or replicons) for the SNP calling.

Note: there must be a space after any comma.
Set to null string to get the largest contig

'''
core_replicon = ""
#core_replicon = "NC_007384, NC_007385"
#core_replicon = "AM412236_4_168118-212711"
#core_replicon = "AM412236_4_168118-212711, ParatyphiA_AKU1"

'''
Mapping To Use:

For Illumina reads, both BWA sampe/samse and Bowtie2 are now available.
If you don't set the mapping, the default is Bowtie mapping
(depending on whether you have pair-ended or single reads).

For Ion Torrent reads, the default is Bowtie2.


            Illumina
            PE    SE    IT
BWA sampe   Y     N     N
BWA samse   N     Y     N
Bowtie2     Y     Y     Y

Note: You cannot use two mappers at the same time!

'''
#mapping = "bwa"
mapping = "bowtie"

'''
Bowtie2 Mapping:
You can change the preset mapping options used with Bowtie2
with the 'bowtie_map_type' option.

Preset options in --end-to-end mode
--very-fast
--fast
--sensitive
--very-sensitive

Preset options in --local mode
--very-fast-local
--fast-local
--sensitive-local
--very-sensitive-local

default: if bowtie_map_type is set to ""
sensitive-local mapping will be used.

See the Bowtie2 Manual for more information on the presets
http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

'''
#end-to-end mode
#bowtie_map_type = "--very-fast"
#bowtie_map_type = "--fast"
#bowtie_map_type = "--sensitive"
#bowtie_map_type = "--very-sensitive"

#local mode
#bowtie_map_type = "--very-fast-local"
#bowtie_map_type = "--fast-local"
bowtie_map_type = "--sensitive-local"
#bowtie_map_type = "--very-sensitive-local"
'''
The default maximum length for bowtie2 to consider pair-ended reads contiguous
is 2000 - you can change this with bowtie_X_value
'''
bowtie_X_value = 2000

'''
BCFtools SNP calling

consensus caller ["c"] (original) or multiallelic caller ["m"] (new bcftools v1+)

Note: multiallelic caller is used only for calling unique SNPs from the BAM files.
The consensus sequences used to populate the allele table based on these SNPs are
still generated using the original consensus caller - this will be changed if/when
a vcf2fq program is available for multiallelic-generated VCFs
'''
SNPcaller = "c"
#SNPcaller = "m"

'''
You can also "remove" any reads: these will be marked as "failed"
This only works during a "merge run"
eg. replace a set of reads with their qc-ed version
    replaceReads ="read_set_2,read_set_3,read_set_24"

'''
replaceReads = ""
#replaceReads = "pool8_tag1"
#replaceReads = "pool1_tag1,pool1_tag3"

'''
Minimum depth of reads for variant filtering
Default value:
        minimum_depth = 5

'''
minimum_depth = 5
'''
HetsVCF
The pipeline filters out heterozygous SNP calls.
To capture these SNPs in the form of a VCF (one per isolate),
set the following to 'True' (Now default behaviour)
'''
#HetsVCF = False
HetsVCF = True

'''
Values for calling the pass/fail and ingroup/outgroup status of strains
suggested (default) values
        cover_fail = 50
        depth_fail = 10
        mapped_fail = 50
        sd_out = 2

'''
cover_fail = 50
depth_fail = 10
mapped_fail = 50
sd_out = 2

'''
Strand Bias Cutoff value
ie. if ABS(DP4[2]-DP4[3])/(DP4[2]+DP4[3]) < strand_bias_cutoff, include the SNP.
Set to >1 to turn strand bias filtering off.
'''
strand_bias_cutoff = 0.8

'''
To switch on or off the checking of percentage of reads mapped
use the following 'check_reads_mapped'.

By default (i.e. set to "") the pipeline will use the largest replicon.
If set to "off", there will be no check for percentage of reads mapped.

Otherwise, give list of the n replicons to be checked, followed by an 'x'
followed by the ratio of the first n-1 replicons
For a single replicon, just put the replicon.
e.g. check_reads_mapped = "rep_1"
or   check_reads_mapped = "rep_1,rep_2,rep_3,x,0.45,0.3"

i.e. rep1 is 45% of the total genome, rep2 is 30% of the total genome,
and rep3 is 25% of the total genome (by default).

Note: there must be no spaces in the list.

'''
check_reads_mapped = ""
#check_reads_mapped = "off"
#check_reads_mapped = "rep_1"
#check_reads_mapped = "rep_1,rep_2,x,0.6"
#check_reads_mapped = "rep_1,rep_2,rep_3,x,0.45,0.3"

#check_reads_mapped = "CP002555,CP002556,x,0.75"

'''
During allele matrix filtering, you can set the conservation level for missing alleles
this is a ratio between 1.0 (100% conservation - remove all SNPs with even one missing allele call)
and 0.0 (0% conservation - remove no SNPs). By default, the pipeline produces the 95% and
0% conservation matrices, with downstream analysis on the 95% matrix.

By entering a different conservation level (e.g. 0.85), both the 95% and 0% matrices
will still be produced, but so too will the 85% matrix (in this example),
and downsteam analysis carried out on this matrix.

'''
conservation = 0.95

'''
Difference Matrix
The pipeline can produce a difference matrix. Currently this is a pairwise difference count.
To get the difference matrix, set the following to 'True'.
'''
DifferenceMatrix = False
#DifferenceMatrix = True


'''
########################
Rubra pipeline variables (do not delete!):
########################
- logDir: the directory where batch queue scripts, stdout and sterr dumps are stored.
- logFile: the file used to log all jobs that are run.
- style: the default style, one of 'flowchart', 'print', 'run', 'touchfiles'. Can be
overridden by specifying --style on the command line.
- procs: the number of python processes to run simultaneously. This determines the
maximum parallelism of the pipeline. For distributed jobs it also constrains the
maximum total jobs submitted to the queue at any one time.
- verbosity: one of 0 (quiet), 1 (normal), 2 (chatty). Can be overridden by specifying
--verbose on the command line.
- end: the desired tasks to be run. Rubra will also run all tasks which are dependencies
of these tasks. Can be overridden by specifying --end on the command line.
- force: tasks which will be forced to run, regardless of timestamps. Can be overridden
by supplying --force on the command line.
- rebuild: one of 'fromstart','fromend'. Whether to calculate which dependencies will
be rerun by working back from an end task to the latest up-to-date task, or forward
from the earliest out-of-date task. 'fromstart' is the most conservative and
commonly used as it brings all intermediate tasks up to date.

'''
pipeline = {
    "logDir": "log",
    "logFile": "pipeline.log",
    "style": "print",
    "procs": 20,
    "paired": True,
    "verbose": 1,
    "end": ["deleteDir"],
    "force": [],
    "rebuild" : "fromstart"
}
stageDefaults = {
    "distributed": True,
    "walltime": "01:00:00",
    "memInGB": 4,
    "queue": "normal",

    "modules": [
         "python/2.7.15-gcc5",
         "bwa/0.7.17-gcc5",
         "samtools/1.9-gcc5",
         "bcftools/1.8",
         "ea-utils/1.1.2-gcc5",
         "fasttree/2.1.10",
         "bowtie2/2.2.9"
    ]
}

stages = {
    "makeDir": {
        "walltime": "00:10:00",
	"queue": "shortq",
        "command": "python makeDir.py %out %sequence_list"
    },
    "copyRef": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "cp %ref %newRef"
    },
    "makeRef": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "python convertGenbankToFasta.py %gen %ref"
    },
    "buildBowtieIndex": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "bowtie2-build %ref %base"
    },
    "alignBowtiePE": {
        "walltime": "03:00:00",
# large file size (any read set >800MB)
#        "walltime": "08:00:00",
        "command": "bowtie2 %type -x %ref_base -1 %seq1 -2 %seq2 -X %Xvalue | samtools view -ubS - | samtools sort - -o %out"
    },
    "alignBowtie": {
        "walltime": "03:00:00",
# large file size (any read set >800MB)
#        "walltime": "06:00:00",
        "command": "bowtie2 %type -x %ref_base -U %seq | samtools view -ubS - | samtools sort - -o %out"
    },
    "buildBWAIndex": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "bwa index -a is %ref"
    },
    "alignSequence": {
# large file size (any read set >800MB)
#        "walltime": "06:00:00",
        "command": "bwa aln %ref %seq > %out"
    },
    "alignBWAPE": {
        "walltime": "03:00:00",
# large file size (any read set >800MB)
#        "walltime": "06:00:00",
        "command": "bwa sampe %ref %align1 %align2 %seq1 %seq2 | samtools view -ubS - | samtools sort - -o %out"
    },
    "alignBWASE": {
        "walltime": "03:00:00",
# large file size (any read set >800MB)
#        "walltime": "06:00:00",
        "command": "bwa samse %ref %align %seq | samtools view -ubS - | samtools sort - -o %out"
    },
    "checkBam": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "python checkBam.py %type %bam %seq"
    },
    "indexBam": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "samtools index %bam"
    },
    "filterUnmapped": {
        "walltime": "00:15:00",
        "queue": "shortq",
        "command": "samtools view -hub -F 4 %bam | samtools sort - -o %out"
    },
    "getSamStats": {
        "walltime": "00:20:00",
        "queue": "shortq",
        "command": "sam-stats -A -B %bam > %out"
    },
    "indexRef": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "samtools faidx %ref"
    },
    "callRepSNPs": {
        "walltime": "01:00:00",
# large file size (any read set >800MB)
#        "walltime": "06:00:00",
        "command": "samtools mpileup -u -t DP -f %ref %bam -r %replicon | bcftools call -O b %option - > %out"
    },
    "checkpoint": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "python checkpoint.py %outTemp %stage"
    },
    "getConsensus": {
        "walltime": "01:00:00",
# large file size (any read set >800MB)
#        "walltime": "06:00:00",
        "command": "samtools mpileup -q 20 -ugB -f %ref %bam | bcftools call -c - | vcfutils.pl vcf2fq > %output"
    },
    "getCoverage": {
        "walltime": "01:00:00",
# large file size (any read set >800MB)
#        "walltime": "06:00:00",
        "command": "samtools mpileup %bam | cut - -f 1-4 > %out"
    },
    "getCoverByRep": {
        "command": "python getCoverByRep.py %ref %coverage %out"
    },
    "q30VarFilterPE": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "bcftools view -i 'ABS(DP4[2]-DP4[3])/(DP4[2]+DP4[3]) < %bias_cutoff ' %bcfFile | vcfutils.pl varFilter -d %min -D %cover -Q 30 > %out"
    },
    "q30VarFilter": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "bcftools view %bcfFile | vcfutils.pl varFilter -d %min -D %cover -Q 30 > %out"
    },
    "finalFilter": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "python finalFilter.py %vcfFile %out %het %flag"
    },
    "getVcfStats": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "python getVcfStats.py %vcfFile %out"
    },
    "deriveRepStats": {
        "walltime": "00:20:00",
        "queue": "shortq",
        "command": "python deriveRepStats.py %coverFile %replicon %depth %cover %runType %map %check"
    },
    "deriveAllStats": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "python deriveAllStats.py %coverFile"
    },
    "collateRepStats": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "python collateRepStats.py %ref %in %replicon %multiplier %runType %sequence_list"
    },
    "collateAllStats": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "python collateAllStats.py %ref %in %path %sequence_list"
    },
    "mergeOutputs": {
        "command": "cp %inputBam %outDirBam && cp %inputIndex %outDirBam && cp %inputVcf %outDirVcf"
    },
    "mergeAllStats": {
        "command": "python mergeAllStats.py %newStats %mergeDir"
    },
    "mergeRepStats": {
        "command": "python mergeRepStats.py %newStats %multiplier %replace %mergeDir %runType"
    },
    "getRepSNPList": {
        "command": "python getRepSNPList.py %in %replicon %out"
    },
    "deriveAllRepGeneCover": {
       "walltime": "00:15:00",
        "queue": "shortq",
       "command": "python deriveAllRepGeneCover.py %outDir %genbank %in"
    },
    "collateAllRepGeneCover": {
# large data sets (more than 150 samples)
#        "walltime": "03:00:00",
        "command": "python collateAllRepGeneCover.py %inDir %outDir %refName %sequence_list"
    },
    "mergeAllRepGeneCover": {
        "walltime": "00:10:00",
        "queue": "shortq",
# large data sets
#        "walltime": "03:00:00",
#	"queue": "normal",
        "command": "python mergeAllRepGeneCover.py %inDir %outDir %refName %sequence_list"
    },
    "parseGeneContent": {
        "walltime": "00:10:00",
        "queue": "shortq",
# large data sets
#        "walltime": "01:00:00",
#        "queue": "normal",
        "command": "python parseGeneContent.py -g %input -d %input2 -s %out -o %out2"
    },
    "deriveRepAlleleMatrix": {
# large number of SNPs (20000+)
#        "walltime": "02:00:00",
        "command": "python deriveRepAlleleMatrix.py %in %out %ref %replicon %consensus %repStats %merge_prefix"
    },
    "collateRepAlleleMatrix": {
        "command": "python collateRepAlleleMatrix.py %in %out %sequence_list %rep_name"
    },
    "getDifferenceMatrix": {
        "walltime": "00:10:00",
        "queue": "shortq",
# large data sets
#        "walltime": "04:00:00",
#        "queue": "normal",
        "command": "python make_distance_matrix.py -i %in"
    },
    "parseSNPs": {
# large data sets
#        "memInGB": 8,
        "command": "python parseSNPtable.py -m cons,aln,coding -s %input -c %conservation -r %genbank -q %replicon -d %dir"
    },
    "parseSNPsNoGBK": {
# large data sets
#        "memInGB": 8,
        "command": "python parseSNPtable.py -m cons,aln -s %input -c %conservation -d %dir"
    },
    "makeTree": {
        "walltime": "01:00:00",
# large data sets
#        "walltime": "06:00:00",
#        "memInGB": 8,
        "command": "FastTreeDbl -gtr -gamma -nt %input > %output"
   },
    "makeNoTree": {
        "command": "python make_no_tree.py %input %out"
    },
    "deleteDir": {
        "walltime": "00:10:00",
        "queue": "shortq",
        "command": "rm -rf %directory"
    }
}

#end of config file
