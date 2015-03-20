RedDog V1.0beta 200315
====== 
Authors: David Edwards, Bernie Pope, Kat Holt

Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors 
may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Description: 
This program implements a workflow pipeline for short read length
sequencing analysis, including the read mapping task, through to variant
detection, followed by analyses (SNPs only).

It uses Rubra (https://github.com/bjpop/rubra) based on the
Ruffus library.

It supports parallel evaluation of independent pipeline stages,
and can run stages on a cluster environment.

Note: for Illumina paired-end or single reads, or Ion Torrent single reads.
IMPORTANT: See config file/instructions for input options/requirements

current version:
V1.0b       Any fixes from final testing (DE)
                includes handling of "." in replicon name (DE)
            Add licensing information to all scripts (DE)
            one more post analysis script added [for Gubbins recombination analysis] (DE)
            update parseGenContent.py (P/A matrix based on cover and depth) (DE)

previous versions:
V0.1        converted to vcf output via mpileup instead of depreciated pileup (DE)
..
V0.2        tested version V0.1.1 (DE)
V0.2.1      adding statistic reporting and fixed "success" handling (DE)
V0.3        tested version of V0.2.1 (DE)
V0.3.1      added alternative output paths to aid clean up (DE)
            new name: pipe_VariantDiscovery.py (DE)
            added q20 and q30 mpileups and associated stats collection (DE)
            various cleanup of code/naming conventions used in pipeline (DE)
V0.3.2      added alternative path for IT or PE data analysis (DE)
V0.3.2.1    added alternative path for SE data analysis (DE)
            added minimum read of depths option for variant filtering (DE)
V0.3.2.2    update in various tools (bwa, bamtools and tmap) (DE)
            changes only affects config file
V0.3.3      tested version of V0.3.2.2 (DE)
V0.3.4      added options for qc of IT reads (DE)
V0.3.4.1    added size option for qc of IT reads (DE)
            updated statistics reporting - minimum reads (DE)
            version numbers change (1.x to 0.x) (DE)
V0.3.4.2    added pass/fail to stats reporting (DE)
            added outgroup/ingroup to stats reporting (DE)
V0.3.5      tested version of V0.3.4.2 (DE)
V0.3.5.1    update to various tools (BWA and tmap) (DE)
            ouput folders now created within the pipeline (DE)
            added separate folder for success files within the temp folder (DE)
            added two new output folders (bam and vcf) (DE)
V0.3.5.2    tested version of V0.3.4.1 (DE)
V0.3.5.3    changed to pipe_vda including:
                allow for merging of new sets of reads into a prior run (DE)
                inclusion of analysis pipelets (DE)
                    - pipe_VCFAnalysis and pipe_AllGeneCover
                clean up of temp directory and/or output directory (if merging) (DE)
                changed "type" to "readType" (DE)
                slight change to pipeline order (DE)
V0.4        tested version of V0.3.5.3 (DE)
V0.4.0.1    merging of bams from different read sets of same strain (DE)
                (either during "new run" or "merge run")
            fixed bug in gene cover and depth matrices script (DE)
V0.4.0.2    tested version of V0.4.0.1 (DE)
V0.4.0.3    removal of QC from within pipeline (and testing) (DE)
V0.4.0.4    replace filter.awk with python-based filtering of all hets from Q30 vcfs (DE)
            includes counting removed het SNPS and reporting same in stat.tab (and testing) (DE)
V0.4.0.5    inclusion of parseSNPtable script (alignment, SNP consequences) (KH, DE)
            and tree generation (DE)
V0.4.0.6    corrections to many scripts used by pipeline, including allele matrix calling 
                and downstream effects to pipeline (DE)
            allele matrix calling now uses consensus sequences (DE)
            addition of differences of SNPs as distance matrix in NEXUS format (DE)
            gene cover and depth matrices no longer contain "fails" (DE)
            addition of parseGeneContent script (KH, DE)
            removal of q20 vcfs (and reporting) - not required (DE)
V0.4.0.7    removed duplicated stages from config file (DE)
V0.4.0.7.1  fix to allow sequences from different folders to be analysed in the same run (DE)
V0.4.0.7.2  fix to deriveStats that let some failed reads pass on depth (DE)
V0.4.1      handling of reference with multiple "chromosomes": pangenome mapping (DE)
                - simplest case: new run (no merging of runs or samples)
                - up to stats collection ('collateRepStats', no post-stats analyses)
            add final '/' to output path(s) if missing (DE)
V0.4.2      handling of reference with multiple "chromosomes": phylogenetic mapping (DE)
                - simplest case: new run (no merging of runs or samples)
                - up to stats collection ('collateRepStats', no post-stats analyses)
            added start-up message (DE)
            changed reference entry from GenBank and Fasta formats to GenBank or Fasta format (DE)
                - fasta reference generated from user GenBank reference
            added pre-run checks including 
                - pairs of reads exist before starting PE analysis (DE)
                - check for 'sequence' option - bad pattern entry (DE)
                - valid run and read types are entered (DE)
            zeroing of SAM files when no longer needed (DE)
V0.4.3      added 'post-stats' analyses - pangenome and phylogeny - no genbank (DE)
            added pre-run reporting and run start confirmation (DE)
V0.4.4      added 'post-stats' analyses - pangenome and phylogeny - with genbank (DE)
V0.4.4.1    various small fixes (DE)
V0.4.4.2    more various small fixes (DE)
V0.4.4.3    increased speed of deriveAllRepGeneCover and getCoverByRep (DE)
V0.4.4.4    conversion of pipeline to use SLURMed Rubra (DE)
V0.4.4.5    fix for bug in BWA sampe/samse v0.7.5 (DE)
V0.4.5      renamed pipeline (DE)
            add merging of runs for pangenome and phylogenetic mapping (DE)
            remove single replicon run (DE)
            added bowtie2 to mapping options (all read types) (DE)
            removal of tmap (DE)
            removal and replacement of bamtools (pileup for coverage instead) (DE)
            cleanup of pipeline scripting (amalgamation of repeated stages) (DE)
            converted emboss call to a biopython script (DE)
            add 'check_reads_mapped' variable for multiple replicon runs (DE)
V0.4.5.1    fix for replicon statistics generation for pangenome runs (DE)
V0.4.5.1.1  fix for all statistics generation when no reads map (DE)
V0.4.5.2    check that replicons all have unique names (DE)
            check that output and out_merge_target folders are different (DE)
            check that output folder is not empty string (DE)
            splitting of getRepAlleleMatrix to improve performance (DE)
                includes sequence list generation (start of post-run report)
V0.4.6      update to newer version of parseSNPtable.py (DE)
                - generation of variable and conserved SNP tables (DE)
                - includes of additional option of setting conservation level (DE)
            further early checks that include:
                - name of reference/replicons/isolates won't confuse post-NEXUS analysis (i.e. no '+')
            fix for when a replicon consensus fasta is missing (DE)
                includes new 'warning' file (DE)
            change behaviour of outgroups - reported (also in outgroups.txt fle) but not removed (DE)
            Editing and reorder of options in config file (DE)
V0.4.7      changed 'sequence_list.txt' generation to function (DE)
            added stage counts and check before firing last stage deleteDir (DE)
            added check for isolates/reads with identical names (DE)
            fixes for errors in mergeRepStats and parseSNPtable (DE)
                - latter includes fixes in script to improve performance (DE)
V0.4.8      include post-run report file function (DE)
            test for 'output' folder prior to run (DE)
            default conservation changed to 0.95 (DE)
            consolidated chrom_info functions into pipe_utils (DE)
            further replicon name checking (DE)
            fix for pipe-generated gene 'tags' when missing (DE)
            write cns warning files to outMerge on merge run (DE)
V0.4.9      splitting location of intermediate files in temp folder to improve stability for large runs (DE)
            changed -X switch in bowtie2 PE mapping from 500 (default) to 1500 (DE)
            removal of some redundant scripts and stages in config file (DE)
V0.5        added check for deletion of previous run success file on merge run (DE)
            added checkpoints for better pipeline running - will halt on errors as expected (DE)
            includes changes to complex stages - flagFiles behaviour (DE)
V0.5.1      added -X option for bowtie2 mapping (DE)
            fixed parseGeneContent output (DE)
            updated to fixed and extended parseSNPtable (DE)
            added scripts for tutorial (filterCoords.py, get_cover.py, getRecomb.R, plotTree.R) (DE)
            changed getDifferenceMatrix to optional output in pipe, changed script to take options (DE)
            added option for VCF output of filtered hets (DE)
            implemented changes to run report from user feedback (DE)
V0.5.2      run report provides settings for merge runs (continuity checks) (DE)
                includes more robust 'check_reads_mapped' (DE)
            update to use SAMtools v1+ (DE)
                includes (limited) addition of multiallelic SNP calling option - bcftools (DE)
            added checkpoint_getSamStats to capture failure during initial BAM construction (DE)
            changed bams from glob call to list call (DE)
            various small fixes to some default values (DE)
            fixed getRecomb.R (DE)

next planned updates
V1.0        Offical Public Release of RedDog (on acceptance of paper) 

(Post-release)
    changes to input sequences pattern recognition
    extend parseGeneContent functions 
    add merging of bams for both pangenome and phylogenetic mapping
    reanalysis without mapping
        with/without a GenBank file,
        restore of read sets removed by user,
        merging of prior runs,
        recalculated/user-edited 'stats.txt' option 
        merging of bams
    further analysis options 

NOTE: with a workaround, some reanalysis without mapping IS possible. Email me for details (DE)

If you wish to see other options added, email me (DE) with suggestions:
(I'm not making any promises...)
    davidje at student dot unimelb dot edu dot au

(DE) My (ongoing) thanks to the "alpha-testers" for their feedback and patience.

