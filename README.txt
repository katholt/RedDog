RedDog V1beta.10.4 260719 ("StopBreakingItStephen")
====== 
Copyright (c) 2016 David Edwards, Bernie Pope, Kat Holt
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

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
V1beta.10.4 ("StopBreakingItStephen") Various fixes including stopping mapping to 
            references with ambiguous base calls
            
previous versions:
V1beta.10.3.1 (Calico Cat)
            fix for parseSNPtabe
V1beta.10.3 (Calico Cat)
            fix for SE and IT reads (no stand bias calling atm - need to validate b4 adding)
            update to many programs used (local Helix installation)
            User Manual update
V1beta.10.2 minor corrections for local install and R scripts
V1beta.10.1 fix for minor error in q30VarFilter
V1beta.10   added strand bias test
V0.1        converted to vcf output via mpileup instead of depreciated pileup 
..
V0.2        tested version V0.1.1 
V0.2.1      adding statistic reporting and fixed "success" handling 
V0.3        tested version of V0.2.1 
V0.3.1      added alternative output paths to aid clean up 
            new name: pipe_VariantDiscovery.py 
            added q20 and q30 mpileups and associated stats collection 
            various cleanup of code/naming conventions used in pipeline 
V0.3.2      added alternative path for IT or PE data analysis 
V0.3.2.1    added alternative path for SE data analysis 
            added minimum read of depths option for variant filtering 
V0.3.2.2    update in various tools (bwa, bamtools and tmap) 
            changes only affects config file
V0.3.3      tested version of V0.3.2.2 
V0.3.4      added options for qc of IT reads 
V0.3.4.1    added size option for qc of IT reads 
            updated statistics reporting - minimum reads 
            version numbers change (1.x to 0.x) 
V0.3.4.2    added pass/fail to stats reporting 
            added outgroup/ingroup to stats reporting 
V0.3.5      tested version of V0.3.4.2 
V0.3.5.1    update to various tools (BWA and tmap) 
            ouput folders now created within the pipeline 
            added separate folder for success files within the temp folder 
            added two new output folders (bam and vcf) 
V0.3.5.2    tested version of V0.3.4.1 
V0.3.5.3    changed to pipe_vda including:
                allow for merging of new sets of reads into a prior run 
                inclusion of analysis pipelets 
                    - pipe_VCFAnalysis and pipe_AllGeneCover
                clean up of temp directory and/or output directory (if merging) 
                changed "type" to "readType" 
                slight change to pipeline order 
V0.4        tested version of V0.3.5.3 
V0.4.0.1    merging of bams from different read sets of same strain 
                (either during "new run" or "merge run")
            fixed bug in gene cover and depth matrices script 
V0.4.0.2    tested version of V0.4.0.1 
V0.4.0.3    removal of QC from within pipeline (and testing) 
V0.4.0.4    replace filter.awk with python-based filtering of all hets from Q30 vcfs 
            includes counting removed het SNPS and reporting same in stat.tab (and testing) 
V0.4.0.5    inclusion of parseSNPtable script (alignment, SNP consequences) (KH, DE)
            and tree generation 
V0.4.0.6    corrections to many scripts used by pipeline, including allele matrix calling 
                and downstream effects to pipeline 
            allele matrix calling now uses consensus sequences 
            addition of differences of SNPs as distance matrix in NEXUS format 
            gene cover and depth matrices no longer contain "fails" 
            addition of parseGeneContent script (KH, DE)
            removal of q20 vcfs (and reporting) - not required 
V0.4.0.7    removed duplicated stages from config file 
V0.4.0.7.1  fix to allow sequences from different folders to be analysed in the same run 
V0.4.0.7.2  fix to deriveStats that let some failed reads pass on depth 
V0.4.1      handling of reference with multiple "chromosomes": pangenome mapping 
                - simplest case: new run (no merging of runs or samples)
                - up to stats collection ('collateRepStats', no post-stats analyses)
            add final '/' to output path(s) if missing 
V0.4.2      handling of reference with multiple "chromosomes": phylogenetic mapping 
                - simplest case: new run (no merging of runs or samples)
                - up to stats collection ('collateRepStats', no post-stats analyses)
            added start-up message 
            changed reference entry from GenBank and Fasta formats to GenBank or Fasta format 
                - fasta reference generated from user GenBank reference
            added pre-run checks including 
                - pairs of reads exist before starting PE analysis 
                - check for 'sequence' option - bad pattern entry 
                - valid run and read types are entered 
            zeroing of SAM files when no longer needed 
V0.4.3      added 'post-stats' analyses - pangenome and phylogeny - no genbank 
            added pre-run reporting and run start confirmation 
V0.4.4      added 'post-stats' analyses - pangenome and phylogeny - with genbank 
V0.4.4.1    various small fixes 
V0.4.4.2    more various small fixes 
V0.4.4.3    increased speed of deriveAllRepGeneCover and getCoverByRep 
V0.4.4.4    conversion of pipeline to use SLURMed Rubra 
V0.4.4.5    fix for bug in BWA sampe/samse v0.7.5 
V0.4.5      renamed pipeline 
            add merging of runs for pangenome and phylogenetic mapping 
            remove single replicon run 
            added bowtie2 to mapping options (all read types) 
            removal of tmap 
            removal and replacement of bamtools (pileup for coverage instead) 
            cleanup of pipeline scripting (amalgamation of repeated stages) 
            converted emboss call to a biopython script 
            add 'check_reads_mapped' variable for multiple replicon runs 
V0.4.5.1    fix for replicon statistics generation for pangenome runs 
V0.4.5.1.1  fix for all statistics generation when no reads map 
V0.4.5.2    check that replicons all have unique names 
            check that output and out_merge_target folders are different 
            check that output folder is not empty string 
            splitting of getRepAlleleMatrix to improve performance 
                includes sequence list generation (start of post-run report)
V0.4.6      update to newer version of parseSNPtable.py 
                - generation of variable and conserved SNP tables 
                - includes of additional option of setting conservation level 
            further early checks that include:
                - name of reference/replicons/isolates won't confuse post-NEXUS analysis (i.e. no '+')
            fix for when a replicon consensus fasta is missing 
                includes new 'warning' file 
            change behaviour of outgroups - reported (also in outgroups.txt fle) but not removed 
            Editing and reorder of options in config file 
V0.4.7      changed 'sequence_list.txt' generation to function 
            added stage counts and check before firing last stage deleteDir 
            added check for isolates/reads with identical names 
            fixes for errors in mergeRepStats and parseSNPtable 
                - latter includes fixes in script to improve performance 
V0.4.8      include post-run report file function 
            test for 'output' folder prior to run 
            default conservation changed to 0.95 
            consolidated chrom_info functions into pipe_utils 
            further replicon name checking 
            fix for pipe-generated gene 'tags' when missing 
            write cns warning files to outMerge on merge run 
V0.4.9      splitting location of intermediate files in temp folder to improve stability for large runs 
            changed -X switch in bowtie2 PE mapping from 500 (default) to 1500 
            removal of some redundant scripts and stages in config file 
V0.5        added check for deletion of previous run success file on merge run 
            added checkpoints for better pipeline running - will halt on errors as expected 
            includes changes to complex stages - flagFiles behaviour 
V0.5.1      added -X option for bowtie2 mapping 
            fixed parseGeneContent output 
            updated to fixed and extended parseSNPtable 
            added scripts for tutorial (filterCoords.py, get_cover.py, getRecomb.R, plotTree.R) 
            changed getDifferenceMatrix to optional output in pipe, changed script to take options 
            added option for VCF output of filtered hets 
            implemented changes to run report from user feedback 
V0.5.2      run report provides settings for merge runs (continuity checks) 
                includes more robust 'check_reads_mapped' 
            update to use SAMtools v1+ 
                includes (limited) addition of multiallelic SNP calling option - bcftools 
            added checkpoint_getSamStats to capture failure during initial BAM construction 
            changed bams from glob call to list call 
            various small fixes to some default values 
            fixed getRecomb.R 
V1.0b       Any fixes from final testing 
                includes handling of "." in replicon name 
            Add licensing information to all scripts 
            one more post analysis script added [for Gubbins recombination analysis] 
            update parseGeneContent.py (P/A matrix based on cover and depth) 
V1beta.1    fix for parseSNPTable - reported position of snp in non-coding feature 
V1beta.2    fix for parseSNPTable - improved speed of reading and parsing snp table 
V1beta.3    changed fasttree to raxml 
            added option to stop tree generation,
            or force tree if > 200 isolates 
V1beta.4    added simple check for correct BAM generation 
            added checkpoint for consensus calling 
V1beta.5    added further filter of SNPS in finalFilter 
V1beta.6    changed back to FastTree - precision error in RAxML -m ASC_GTRCAT 
            changed maximum isolates for tree to 500 
            changed checkBam to pass BAMs from simulated reads 
V1beta.7 (BlackCat)
            fixed bug in quality filtering of variant calls 
V1beta.8    tutorial update 
V1beta.9    local system update 
            manual update 

