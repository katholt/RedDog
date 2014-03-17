'''
Configuration file for RedDog.py V0.4.5.2
-------------------------------

Reference and sequences (from VR0082 shared directory)

Reference can be GenBank or fasta format - if GenBank format, this will be converted
to a fasta version for mapping.
 
If the GenBank file is not given, the gene cover and depth matrices for genes
will not be generated and nor will SNP consequences. 

If you don't have the GenBank record, or don't want the above matrices 
to be generated, enter a fasta format reference instead.

'''
#Test Sets
#reference = "/vlsci/VR0082/shared/pipeline_test_data/reference/NC_007384.gbk"
reference = "/vlsci/VR0082/shared/pipeline_test_sets/reference/NC_007384_with_plasmid.gbk"
sequences = "/vlsci/VR0082/shared/pipeline_test_sets/illumina/shigella/*.fastq.gz"
#sequences = "/vlsci/VR0082/shared/pipeline_test_sets/illumina/shigella/extra/*.fastq.gz"

# You can now also combine sequences from different folders into the same run...
#sequences = ["/vlsci/VR0082/shared/pipeline_test_data/illumina_pe/*.fastq.gz", "/vlsci/VR0082/shared/pipeline_test_data/illumina_pe/extra/*.fastq.gz"]

#reference = "/vlsci/VR0082/shared/pipeline_test_sets/reference/DT104.fasta"
#sequences = "/vlsci/VR0082/shared/pipeline_test_data/salmonella/*in.iontor.fastq.gz"


'''
Choose type of input sequences: IT for ion torrent (single reads), 
                                PE for Illumina pair-end reads
				                        or SE for Illumina single-end reads

readType = "PE" or "SE" or "IT"

Each read type has a particular pattern you need to follow for use in the pipeline:

SE: *.fastq.gz                      (i.e. must be gzipped)
PE: *_1.fastq.gz and *_2.fastq.gz   (i.e. must be gzipped with forward and reverse in separate files)
IT: *_in.iontor.fastq.gz            (i.e. must be gzipped!!!!)

'''
readType = "PE"
#readType = "SE"
#readType = "IT"

'''
Run Type - if set to a null string, the number of replicons in the reference will determine
the run type:
    1 - 100 replicons (e.g. reference genome + plasmids + phage)   - phylogeny run type
    > 100 replicons (e.g. multifasta pangenome)                    - pangenome run type

The user can override the run type, by setting it below. If there is only one replicon,
the pipeline will default to a standard (single reference) run. The run types are described in 
more detail in the instructions.

'''
#runType = ""
runType = "pangenome"
#runType = "phylogeny"

'''
For a pangenome run, the SNPs will only be called for the largest replicon - this assumes
the core genome is in this replicon. The user can define an alternative replicon
(or replicons) for the SNP calling.

Set to null string to get the largest contig, or for any other run type
'''
core_replicon = ""
#core_replicon = "AM412236_4_168118-212711"
#core_replicon = "AM412236_4_168118-212711, ParatyphiA_AKU1"
#Salmonella Typhimurium STm135
#core_replicon = "NC_016810, NC_017720"

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
Output directory:
full path name including final "/"

VR0082 users: Make sure this is to a directory in the shared folder!!!
'''
output = "/vlsci/VR0082/shared/davide/pipe_test_out/mapping/NC_007384_pan/"

'''
Directory to merge output with (out_merge_target):
 
When running new analysis set to null string.
 
Otherwise set to the directory you want to merge with.
 
You can only merge a prior run with a new run (not two prior runs)
This merge target folder must have the bams and indexes in one sub-folder (/bam)
and the vcfs in another (/vcf). 
There also must be a sequence_list.txt file - i.e. V0.4.5.2+ format.
 
The 'output' folder (see above) for a merge run should NOT exist prior to the run,
and will be deleted at completion of the pipeline.

Set to empty string for no merging (i.e. new run).

Note: the pipeline can no longer merge 'single' run types (those that use 'stats.tab).
If you really need to do so, make use of v0.4.4.4 of the pipeline... 
'''
out_merge_target = ""
#out_merge_target = "/vlsci/VR0082/shared/davide/pipe_test_out/mapping/NC_007384_phy/"

'''
You can also "replace" any reads: these will be marked as "failed"
This only works during a "merge run"
eg. replace a set of reads with their qc-ed version
    replaceReads ="'read_set_2', 'read_set_3'"
'''
replaceReads = ""
#replaceReads = "'pool8_tag1', 'pool8_tag2'"

'''
Minimum depth of reads for variant filtering
Default value:
        minimum_depth = 5
'''
minimum_depth = 5

'''
Values for calling the pass/fail and ingroup/outgroup status of strains
suggested (default) values for standard run (single reference genome)
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
To switch on or off the checking of percentage of reads mapped
use the following 'check_reads_mapped'.

By default (i.e. set to "") the pipeline will use the largest replicon.
If set to "off", there will be no check for percentage of reads mapped.

Otherwise, give list of the n replicons to be checked, followed by an 'x'
followed by the ratio of the first n-1 replicons
For a single replicon, just put the replicon.
e.g. 
check_reads_mapped = "rep_1"

or

check_reads_mapped = "rep_1,rep_2,rep_3,x,0.45,0.3"
i.e. rep1 is 45% of the total genome, rep2 is 30% of the total genome,
and rep3 is 25% of the total genome (by default). 

Note: there must be no spaces in the list... 
'''
check_reads_mapped = ""
#check_reads_mapped = "off"
#check_reads_mapped = "rep_1"
#check_reads_mapped = "rep_1,rep_2,x,0.6"
#check_reads_mapped = "rep_1,rep_2,rep_3,x,0.45,0.3"

#check_reads_mapped = "CP002555,CP002556,x,0.75"

#To be added: an option for the user defining the outgroups in a set
#eg. outgroups = "pool1_tag4", "pool10_tag6"
#outgroups = ""
# Not Yet Implemented - don't set!

'''
Rubra variables:
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
    "procs": 100,
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
    "queue": None,
    "modules": [
        "python-gcc/2.7.5",
        "bwa-intel/0.6.2",
        "samtools-intel/0.1.19",
        "eautils-gcc/1.1.2",
        "fasttree-intel/2.1.7",
        "bowtie2-intel/2.1.0"
    ]
}
stages = {
    "makeDir": {
        "walltime": "00:10:00",
        "command": "mkdir -p %dir1 %dir2 %dir3"
    },
    "copyRef": {
        "walltime": "00:10:00",
        "command": "cp %ref %newRef"
    },
    "makeRef": {
        "walltime": "00:10:00",
        "command": "python convertGenbankToFasta.py %gen %ref"
    },
    "buildBowtieIndex": { 
        "walltime": "00:10:00",
        "command": "bowtie2-build %ref %base"
    },
    "alignBowtiePE": {
        "walltime": "03:00:00",
        "command": "bowtie2 %type -x %ref_base -1 %seq1 -2 %seq2 | samtools view -ubS - | samtools sort - %out"
    },
    "alignBowtie": {
        "walltime": "03:00:00",
        "command": "bowtie2 %type -x %ref_base -U %seq | samtools view -ubS - | samtools sort - %out"
    },
    "buildBWAIndex": {
        "walltime": "00:10:00",
        "command": "bwa index -a is %ref"
    },
    "alignSequence": {
        "command": "bwa aln %ref %seq > %out"
    },
    "alignBWAPE": {
        "walltime": "03:00:00",
        "command": "bwa sampe %ref %align1 %align2 %seq1 %seq2 | samtools view -ubS - | samtools sort - %out"
    },
    "alignBWASE": {
        "walltime": "02:00:00",
        "command": "bwa samse %ref %align %seq | samtools view -ubS - | samtools sort - %out"
    },
    "indexBam": {
        "walltime": "00:10:00",
        "command": "samtools index %bam"
    },
    "filterUnmapped": {
        "walltime": "00:15:00",
        "command": "samtools view -hub -F 4 %bam | samtools sort - %out"
    },
    "getSamStats": {
        "walltime": "00:20:00",
        "command": "sam-stats -A -B %bam > %out"
    },
    "indexRef": {
        "walltime": "00:10:00",
        "command": "samtools faidx %ref"
    },
    "callSNPs": {
        "walltime": "03:00:00",
        "command": "samtools mpileup -uD -f %ref %bam  | bcftools view -bvcg - > %out"
    },
    "callRepSNPs": {
        "walltime": "01:00:00",
        "command": "samtools mpileup -uD -f %ref %bam -r %replicon | bcftools view -bvcg - > %out"
    },
    "getConsensus": {
        "walltime": "01:00:00",
        "command": "samtools mpileup -q 20 -uB -f %ref %bam | bcftools view -c - | vcfutils.pl vcf2fq > %output"
    },
    "getCoverage": {
        "walltime": "01:00:00",
        "command": "samtools mpileup %bam | cut - -f 1-4 > %out"
    },
    "averageCoverage": {
        "command": "python averageCoverage.py %coverage %minDepth %out"
    },
    "getCoverByRep": {
        "command": "python getCoverByRep.py %ref %coverage %out"
    },
    "q30VarFilter": {
        "walltime": "00:10:00",
        "command": "bcftools view %bcfFile | vcfutils.pl varFilter -d %min -D %cover -Q 30 > %out"
    },
    "finalFilter": {
        "walltime": "00:10:00",
        "command": "python finalFilter.py %vcfFile %out"
    },
    "getVcfStats": {
        "walltime": "00:10:00",
        "command": "python getVcfStats.py %vcfFile %out"
    },
    "deriveRepStats": {
        "walltime": "00:10:00",
        "command": "python deriveRepStats.py %coverFile %replicon %depth %cover %runType %map %check"
    },
    "deriveAllStats": {
        "walltime": "00:10:00",
        "command": "python deriveAllStats.py %coverFile"
    },
    "collateRepStats": {
        "walltime": "00:10:00",
        "command": "python collateRepStats.py %ref %in %replicon %multiplier %runType"
    },
    "collateAllStats": {
        "walltime": "00:10:00",
        "command": "python collateAllStats.py %ref %in"
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
       "command": "python deriveAllRepGeneCover.py %outDir %genbank %in"
    },
    "collateAllRepGeneCover": {
        "walltime": "00:10:00",
# large data sets
#        "walltime": "03:00:00",
        "command": "python collateAllRepGeneCover.py %inDir %outDir %refName"
    },
    "mergeAllRepGeneCover": {
        "walltime": "00:10:00",
# large data sets
#        "walltime": "03:00:00",
        "command": "python mergeAllRepGeneCover.py %inDir %outDir %refName"
    },
    "parseGeneContent": {
        "walltime": "00:10:00",
# large data sets
#        "walltime": "01:00:00",
        "command": "python parseGeneContent.py -g %input -o %out -s %out2"
    },
#replacing this step
    "getRepAlleleMatrix": {
# large data sets
#        "walltime": "02:00:00",
#        "memInGB": 8,
        "command": "python getRepAlleleMatrix.py %in %out %ref %replicon"
    },
#with these next two
    "deriveRepAlleleMatrix": {
# large data sets
#        "walltime": "02:00:00",
#        "memInGB": 8,
        "command": "python deriveRepAlleleMatrix.py %in %out %ref %replicon %consensus %repStats"
    },
    "collateRepAlleleMatrix": {
# large data sets
#        "walltime": "02:00:00",
#        "memInGB": 8,
        "command": "python collateRepAlleleMatrix.py %in %out %length"
    },
    "getDifferenceMatrix": {
        "walltime": "00:10:00",
# large data sets
#        "walltime": "12:00:00",
        "command": "python make_distance_matrix.py %in"
    },
    "parseSNPs": {
# large data sets
#        "walltime": "3:00:00:00",
#        "memInGB": 16,
        "command": "wDir=\\\"`pwd`\\\" && cd %dir && python $wDir/parseSNPtable.py -m aln,coding -r %genbank -s %input"
    },
    "parseSNPsNoGBK": {
        "walltime": "00:10:00",
# large data sets
#        "walltime": "08:00:00",
#        "memInGB": 64,
        "command": "wDir=\\\"`pwd`\\\" && cd %dir && python $wDir/parseSNPtable.py -m aln -s %input"
    },
    "makeTree": {
        "walltime": "00:15:00",
# large data sets
#        "walltime": "18:00:00",
#        "memInGB": 16,
        "command": "FastTree -gtr -gamma -nt %input > %output"
    },
    "deleteDir": {
        "walltime": "00:10:00",
        "command": "rm -rf %directory"
    }
}
