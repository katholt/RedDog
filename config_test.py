'''
Configuration file for RedDog.py V0.5.2
-------------------------------
Essential pipeline variables.
'''
reference = "/vlsci/VR0082/shared/pipeline_test_sets/reference/NC_007384_with_plasmid.gbk"

#sequences = "/vlsci/VR0082/shared/pipeline_test_sets/illumina/shigella/*.fastq.gz"
sequences = "/vlsci/VR0082/shared/pipeline_test_sets/illumina/shigella/extra/*.fastq.gz"

#output = "/scratch/VR0082/workspace/mapping/v052_test"
output = "/scratch/VR0082/workspace/mapping/v052_test_merge"

#out_merge_target = ""
out_merge_target = "/scratch/VR0082/workspace/mapping/v052_test"


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
        "samtools-intel/1.1",
        "bcftools-intel/1.1",
        "eautils-gcc/1.1.2",
        "fasttree-intel/2.1.7",
        "bowtie2-intel/2.2.3"
    ]
}
stages = {
    "makeDir": {
        "walltime": "00:10:00",
        "command": "python makeDir.py %out %sequence_list"
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
# large file size (any read set >800MB)
#        "walltime": "06:00:00",
        "command": "bowtie2 %type -x %ref_base -1 %seq1 -2 %seq2 -X %Xvalue | samtools view -ubS - | samtools sort - %out"
    },
    "alignBowtie": {
        "walltime": "03:00:00",
# large file size (any read set >800MB)
#        "walltime": "06:00:00",
        "command": "bowtie2 %type -x %ref_base -U %seq | samtools view -ubS - | samtools sort - %out"
    },
    "buildBWAIndex": {
        "walltime": "00:10:00",
        "command": "bwa index -a is %ref"
    },
    "alignSequence": {
# large file size (any read set >800MB)
#        "walltime": "02:00:00",
        "command": "bwa aln %ref %seq > %out"
    },
    "alignBWAPE": {
        "walltime": "03:00:00",
# large file size (any read set >800MB)
#        "walltime": "06:00:00",
        "command": "bwa sampe %ref %align1 %align2 %seq1 %seq2 | samtools view -ubS - | samtools sort - %out"
    },
    "alignBWASE": {
        "walltime": "03:00:00",
# large file size (any read set >800MB)
#        "walltime": "06:00:00",
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
    "callRepSNPs": {
        "walltime": "01:00:00",
# large file size (any read set >800MB)
#        "walltime": "03:00:00",
        "command": "samtools mpileup -u -t DP -f %ref %bam -r %replicon | bcftools call -O b %option - > %out"
    },
    "checkpoint": {
        "walltime": "00:10:00",
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
        "command": "time samtools mpileup %bam | cut - -f 1-4 > %out"
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
        "command": "python finalFilter.py %vcfFile %out %het %flag"
    },
    "getVcfStats": {
        "walltime": "00:10:00",
        "command": "python getVcfStats.py %vcfFile %out"
    },
    "deriveRepStats": {
        "walltime": "00:20:00",
        "command": "python deriveRepStats.py %coverFile %replicon %depth %cover %runType %map %check"
    },
    "deriveAllStats": {
        "walltime": "00:10:00",
        "command": "python deriveAllStats.py %coverFile"
    },
    "collateRepStats": {
        "walltime": "00:10:00",
        "command": "python collateRepStats.py %ref %in %replicon %multiplier %runType %sequence_list"
    },
    "collateAllStats": {
        "walltime": "00:10:00",
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
       "command": "python deriveAllRepGeneCover.py %outDir %genbank %in"
    },
    "collateAllRepGeneCover": {
        "walltime": "00:10:00",
# large data sets (more than 150 samples)
#        "walltime": "03:00:00",
        "command": "python collateAllRepGeneCover.py %inDir %outDir %refName %sequence_list"
    },
    "mergeAllRepGeneCover": {
        "walltime": "00:10:00",
# large data sets
#        "walltime": "03:00:00",
        "command": "python mergeAllRepGeneCover.py %inDir %outDir %refName %sequence_list"
    },
    "parseGeneContent": {
        "walltime": "00:10:00",
# large data sets
#        "walltime": "01:00:00",
        "command": "python parseGeneContent.py -g %input -s %out -o %out2"
    },
    "deriveRepAlleleMatrix": {
        "command": "python deriveRepAlleleMatrix.py %in %out %ref %replicon %consensus %repStats %merge_prefix"
    },
    "collateRepAlleleMatrix": {
        "command": "python collateRepAlleleMatrix.py %in %out %sequence_list %rep_name"
    },
    "getDifferenceMatrix": {
        "walltime": "00:10:00",
# large data sets
#        "walltime": "04:00:00",
        "command": "python make_distance_matrix.py -i %in"
    },
    "parseSNPs": {
# large data sets
#        "walltime": "03:00:00",
#        "memInGB": 8,
        "command": "python parseSNPtable.py -m cons,aln,coding -s %input -c %conservation -r %genbank -q %replicon -d %dir"
    },
    "parseSNPsNoGBK": {
        "walltime": "00:10:00",
# large data sets
#        "walltime": "03:00:00",
#        "memInGB": 8,
        "command": "python parseSNPtable.py -m cons,aln -s %input -c %conservation -d %dir"
    },
    "makeTree": {
        "walltime": "00:15:00",
# large data sets
#        "walltime": "06:00:00",
#        "memInGB": 8,
        "command": "FastTree -gtr -gamma -nt %input > %output"
    },
    "deleteDir": {
        "walltime": "00:10:00",
        "command": "rm -rf %directory"
    }
}

#end of config file
