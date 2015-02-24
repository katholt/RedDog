#!/bin/env python

'''
RedDog V0.5.2 240215
====== 
Authors: David Edwards, Bernie Pope, Kat Holt
License: none as yet...

Description: 

This program implements a worklow pipeline for short read length
sequencing analysis, including the read mapping task, through to variant
detection, followed by analyses (SNPs only).

It uses Rubra (https://github.com/bjpop/rubra) based on the
Ruffus package (v2.2).

It supports parallel evaluation of independent pipeline stages,
and can run stages on a cluster environment.

Note: for Illumina paired-end or single reads, or Ion Torrent single reads.
IMPORTANT: See config file/instructions for input options/requirements

Version History: See ReadMe.txt

'''
from ruffus import *
import os
import shutil
import sys
import glob
from rubra.utils import pipeline_options
from rubra.utils import (runStageCheck, splitPath)
from pipe_utils import (isGenbank, isFasta, chromInfoFasta, chromInfoGenbank, getValue, getCover, make_sequence_list, getSuccessCount, make_run_report, get_run_report_data)

version = "V0.5.2"

modules = pipeline_options.stageDefaults['modules']

# determine the reference file,
# list of sequence files, and list of chromosmes.
try:
    reference = pipeline_options.reference
except:
    print "\nNo Reference supplied"
    print "Pipeline Stopped: please supply a reference in the config file\n"
    sys.exit()    

#check that a reference has been given in the options file
if reference == "":
    print "\nNo Reference supplied"
    print "Pipeline Stopped: please supply a reference in the config file\n"
    sys.exit()    

#check whether reference file exists
if not os.path.exists(reference):
    print "\nReference supplied does not exist"
    print "Pipeline Stopped: please check supplied reference\n"
    sys.exit()    

#check whether reference is in FASTA or GenBank format
if isFasta(reference):
    refGenbank = False
    replicons = chromInfoFasta(reference)

elif isGenbank(reference):
    refGenbank = True
    replicons = chromInfoGenbank(reference)
else:
    print "\nReference not in GenBank or FASTA format"
    print "Pipeline Stopped: please check your reference\n"
    sys.exit()

# check that reference name does not in contain "+"
if reference.find("+") != -1:
    print "\nReference has an illegal character in the name ('+') "
    print "Pipeline Stopped: please change the name of the reference\n"
    sys.exit()

# check that replicons in reference have non-unique names
if len(replicons)>1:
    for i in range(len(replicons)-1):
        for j in range (i+1, len(replicons)):
            if replicons[i][0] != replicons[j][0]:
                pass
            else:
                print "\nReference has replicons with non-unique names: " + replicons[i][0]
                print "Pipeline Stopped: please check your reference\n"
                sys.exit()

# check that replicon names in reference do not in contain ":" or "+"
for i in range(len(replicons)):
    if replicons[i][0].find(":") != -1 or replicons[i][0].find("+") != -1 or replicons[i][0].find("|") != -1:
        print "\nReference has replicon with an illegal character ('|', ':' or '+'): " + replicons[i][0]
        print "Pipeline Stopped: please change the name of the replicon in the reference\n"
        sys.exit()

(refPrefix, refName, refExt) = splitPath(reference)

try:
    sequencePatterns = pipeline_options.sequences
except:
    print "\nNo Sequences option supplied"
    print "Pipeline Stopped: please include 'sequences' in the config file\n"
    sys.exit()    

try:
    runType = pipeline_options.runType
except:
    runType = ""

try:
    core_replicon = pipeline_options.core_replicon
except:
    core_replicon = ""

#check that a 'known' runType has been supplied
if runType != "":
    if runType == "pangenome" or runType == "phylogeny":
        pass
    else:
        print "\nUnrecognised run type"
        print "Pipeline Stopped: please check 'runType' in the config file\n"
        sys.exit()

# if no runType entered, work out default runType on number of replicons in reference
if runType == "":
    if len(replicons) > 100:
        runType = "pangenome"
    elif len(replicons) > 0:
        runType = "phylogeny"
if len(replicons) < 1:
    print "\nNo chromosomes found in reference"
    print "Pipeline Stopped: please check 'reference' in the config file\n"
    sys.exit()

repliconNames = []
longest_replicon_length = 0
for repliconName, repliconLength in replicons:
    repliconNames.append(repliconName)
    if int(repliconLength) > longest_replicon_length:
        longest_replicon_length = int(repliconLength)

core_replicons = []
if runType == "pangenome":
    if core_replicon == "":
        for repliconName, repliconLength in replicons:
            if int(repliconLength) == longest_replicon_length:
                core_replicon = repliconName
    core_replicons = core_replicon.split(', ')

sequences = []
if type(sequencePatterns) == list:
    for pattern in sequencePatterns:
        sequences += glob.glob(pattern)
else:
    sequences = glob.glob(sequencePatterns)

if sequences == []:
    if sequencePatterns == "":
        print "\nNo mapping - analysis only"
        print "Pipeline Stopped: Not Yet Available\n"
        sys.exit()
    else:
        print "\nNo matching sequences found"
        print "Pipeline Stopped: please check 'sequences' in the config file\n"
        sys.exit()

for sequence in sequences:
    if sequence.find(":") != -1 or sequence.find("+") != -1:
        print "\nA read set has an illegal character (':' or '+'): " + sequence
        print "Pipeline Stopped: please change the name of the read set\n"
        sys.exit()

try:
    readType = pipeline_options.readType
except:
    readType = 'PE'    

if readType == 'IT' or readType == 'PE' or readType == 'SE':
    pass
elif readType == '':
    readType = 'PE'
else:
    print "\nUnrecognised read type"
    print "Pipeline Stopped: please check 'readType' in the config file\n"
    sys.exit()

mapping_out = ""
try:
    mapping = pipeline_options.mapping
except:
    mapping = 'bowtie'

if mapping == "":
    mapping = 'bowtie'

if mapping == 'bwa' and readType == 'SE':
    mapping_out = 'BWA V0.6.2 samse'
elif mapping == 'bwa' and readType == 'PE':
    mapping_out = 'BWA V0.6.2 sampe'
elif mapping == 'bowtie':
    mapping_out = 'Bowtie2 V2.2.3'
else:
    print "\nUnrecognised mapping option"
    print "Pipeline Stopped: please check 'mapping' in the config file\n"
    sys.exit()

try:
    SNPcaller = pipeline_options.SNPcaller
    if not(SNPcaller != 'c' or SNPcaller != 'm'):
        print "\nUnrecognised SNPcaller option"
        print "Pipeline Stopped: please check 'SNPcaller' in the config file\n"
    else:
        SNPcaller = '-' + SNPcaller
except:
    SNPcaller = '-c'

sequence_list = []
for sequence in sequences:
    (prefix, name, ext) = splitPath(sequence)
    if readType == "IT":
        sequence_list.append(name[:-16])
    elif readType == "PE":
        if name.find('_1.fastq') != -1 and name[:-8] not in sequence_list:
            sequence_list.append(name[:-8])
        if name.find('_2.fastq') != -1 and name[:-8] not in sequence_list:
            sequence_list.append(name[:-8])            
    else:
        sequence_list.append(name[:-6])

if readType == 'PE':
    missing_pairs = []
    for seqName in sequence_list:
        test_1 = False
        test_2 = False
        for sequence in sequences:
            if sequence.find(seqName+'_1.f') != -1:
                test_1 = True
            if sequence.find(seqName+'_2.f') != -1:
                test_2 = True
        if test_1 == False or test_2 == False:
            missing_pairs.append(seqName)
    if missing_pairs != []:
        print "\nNot all sequence sets have pairs:"
        for seq in missing_pairs:
            print seq
        print "Pipeline Stopped: please fix sequence pairs\n"
        sys.exit()

if mapping == 'bowtie':
    try:
        bowtie_map_type = pipeline_options.bowtie_map_type
    except:
        bowtie_map_type = '--sensitive-local'
    if bowtie_map_type == "":
        bowtie_map_type = '--sensitive-local'
    if (bowtie_map_type == "--very-fast" or
        bowtie_map_type == "--fast" or
        bowtie_map_type == "--sensitive" or
        bowtie_map_type == "--very-sensitive" or
        bowtie_map_type == "--very-fast-local" or
        bowtie_map_type == "--fast-local" or
        bowtie_map_type == "--sensitive-local" or
        bowtie_map_type == "--very-sensitive-local"):
        pass        
    else:
        print "\nUnrecognised Bowtie2 mapping option"
        print "Pipeline Stopped: please check 'bowtie_map_type' in the config file\n"
        sys.exit()
else:
    bowtie_map_type = "-"

try:
    bowtie_X_value = pipeline_options.bowtie_X_value
except:
    bowtie_X_value = 2000

try:
    minDepth = pipeline_options.minimum_depth
except:
    minDepth = 5

try:
    HetsVCF = pipeline_options.HetsVCF
    if HetsVCF != True and HetsVCF != False:
        print "\nUnrecognised HetsVCF option"
        print "Pipeline Stopped: please check 'HetsVCF' in the config file\n"
        sys.exit()
except:
    HetsVCF = False

try:
    coverFail = pipeline_options.cover_fail
except:
    coverFail = 50

try:
    depthFail = pipeline_options.depth_fail
except:
    depthFail = 10

try:
    mappedFail = pipeline_options.mapped_fail
except:
    mappedFail = 50

try:
    sdOutgroupMultiplier = pipeline_options.sd_out
except:
    sdOutgroupMultiplier = 2

try:
    check_reads_mapped = pipeline_options.check_reads_mapped
except:
    check_reads_mapped = ""

if check_reads_mapped == "":
    for repliconName, repliconLength in replicons:
        if int(repliconLength) == longest_replicon_length:
            check_reads_mapped = repliconName
elif check_reads_mapped != 'off':
    if check_reads_mapped.find(',') == -1:
        if check_reads_mapped not in repliconNames:
            print "\nUnrecognised replicon used in 'check_reads_mapped': " + check_reads_mapped
            print "Pipeline Stopped: please check 'check_reads_mapped' in the config file\n"
            sys.exit()
    else:
        found_x = False
        final_ratio = 1.0
        list_of_replicons = []
        check_reps = check_reads_mapped.split(',')
        for item in check_reps:
            if item != 'x':
                if found_x != True:
                    list_of_replicons.append(item)
                else:
                    final_ratio -= float(item)
            else:
                found_x = True
        if final_ratio <= 0:
            print "\nFinal ratio in 'check_reads_mapped' less than or equal to zero"
            print "Pipeline Stopped: please check ratio(s) in 'check_reads_mapped' in the config file\n"
            sys.exit()
        list_of_replicons_done = []
        for replicon in list_of_replicons:
            if replicon not in repliconNames:
                print "\nUnrecognised replicon in 'check_reads_mapped': " + replicon
                print "Pipeline Stopped: please check 'check_reads_mapped' in the config file\n"
                sys.exit()
            elif replicon not in list_of_replicons_done:
                list_of_replicons_done.append(replicon)
            else:
                print "\nReplicon repeated in 'check_reads_mapped': " + replicon
                print "Pipeline Stopped: please check 'check_reads_mapped' in the config file\n"
                sys.exit()

try:
    outPrefix = pipeline_options.output
except:
    outPrefix = ""

if outPrefix == "":
    print "\nNo Output folder given"
    print "Pipeline Stopped: please check 'output' in the config file\n"
    sys.exit()

if outPrefix[-1] != '/':
    outPrefix += '/'

outTempPrefix = outPrefix + 'temp/'
outSuccessPrefix = outTempPrefix + 'success/'
outBamPrefix = outPrefix + 'bam/'
outVcfPrefix = outPrefix + 'vcf/'

if os.path.isdir(outPrefix) and not(os.path.exists(outSuccessPrefix + "dir.makeDir.Success")):
    print "\nOutput folder already exists"
    print "Pipeline Stopped: please change 'output' to a new folder\n"
    sys.exit()    

try:
    outMerge = pipeline_options.out_merge_target
except:
    print "\n'out_merge_target' not set"
    print "Pipeline Stopped: please set 'out_merge_target' in the config file\n"
    sys.exit()    

if outPrefix == outMerge:
    print "\nOutput folder and out_merge_target for run are the same"
    print "Pipeline Stopped: please check 'output' and 'out_merge_target' in the config file\n"
    sys.exit()

continuity_test = False
if outMerge != '':
    if outMerge[-1] != '/':
        outMerge += '/'
    outMergeBam = outMerge + 'bam/'
    outMergeVcf = outMerge + 'vcf/'
    if not(os.path.isdir(outMerge)):
        print "\nMerge target folder not found"
        print "Pipeline Stopped: please change 'out_merge_target' to previous output folder\n"
        sys.exit()
    else:
        if not(os.path.isdir(outMergeBam)):
            print "\nBAM folder not found"
            print "Pipeline Stopped: check 'out_merge_target' has BAM folder\n"
            sys.exit()    
        if not(os.path.isdir(outMergeVcf)):
            print "\nVCF folder not found"
            print "Pipeline Stopped: check 'out_merge_target' has VCF folder\n"
            sys.exit()    
    if os.path.exists(outMerge + 'finish.deleteDir.Success'):
        print "\n'out_merge_target' still has 'finish.deleteDir.Success' file"
        print "Pipeline Stopped: please delete this file before restarting\n"
        sys.exit()
    if os.path.exists(outMerge + refName +'_AllStats.txt') != True:
        print "\n'out_merge_target' has no '"+ refName+"_AllStats.txt' file"
        print "Pipeline Stopped: please check you have the correct 'reference' and/or "
        print "\t\t  'out_merge_target' before restarting\n"
        sys.exit()
    merge_run = True
    if os.path.exists(outMerge + refName + '_run_report.txt'):
        (read_history, 
        run_history, 
        old_ref_name, 
        old_ref_format, 
        old_replicon_count, 
        old_core_replicon_list, 
        old_run_type, 
        old_bowtie_preset, 
        old_bowtie_X, 
        old_user_failed_list, 
        old_min_depth, 
        old_cover_fail, 
        old_depth_fail, 
        old_mapped_fail, 
        old_replicon_test_list, 
        old_replicon_percent_list, 
        old_conservation,
        old_sd_out,
        old_replicon_list) = get_run_report_data(outMerge + refName + '_run_report.txt')
        continuity_test = True
    else:
        print "\nMerge Run: No prior run report found"
        print "No continuity tests will be done...\n"
        run_history = '-'
        read_history = '_'
else:
    merge_run = False
    run_history = '-'
    read_history = '-'

try:
    if pipeline_options.replaceReads != "":
        replaceReads = pipeline_options.replaceReads
    else:
        replaceReads = "-"
except:
    replaceReads = "-"

try:
    conservation = float(pipeline_options.conservation)
except:
    conservation = 0.95

if conservation > 1.0 or conservation < 0.0:
    print "\n'conservation' set to value outside parameters"
    print "Pipeline Stopped: please change 'conservation' to a value between 0 and 1\n"
    sys.exit()

try:
    DifferenceMatrix = pipeline_options.DifferenceMatrix
    if DifferenceMatrix != True and DifferenceMatrix != False:
        print "\nUnrecognised DifferenceMatrix option"
        print "Pipeline Stopped: please check 'DifferenceMatrix' in the config file\n"
        sys.exit()
except:
    DifferenceMatrix = False

full_sequence_list = []
if outMerge == '':
    for item in sequence_list:
        full_sequence_list.append(item)
else:
    for item in sequence_list:
        full_sequence_list.append(item)
    try:
        sequence_list_file = open((outMerge + 'sequence_list.txt') , "r")
    except:
        print "\nNo sequence list found"
        print "Pipeline Stopped: please generate new 'sequence_list.txt' file\n"
        sys.exit()
    for line in sequence_list_file:
        full_sequence_list.append(line[:-1])
    sequence_list_file.close()

duplicate_isolate_name = []
for i in range(len(full_sequence_list)-1):
    for j in range(i+1, len(full_sequence_list)):
        if full_sequence_list[i] == full_sequence_list[j]:
            duplicate_isolate_name.append(i)

if duplicate_isolate_name != []:
    print "\nDuplicate reads (isolate) names found"
    print "Pipeline Stopped: please remove/rename duplicates"
    for i in duplicate_isolate_name:
        print full_sequence_list[i]
    print "\n" 
    sys.exit()

full_sequence_list_string = ""
for sequence in full_sequence_list:
    full_sequence_list_string += sequence + ","
full_sequence_list_string.rstrip()

sequence_list_string = ""
for sequence in sequence_list:
    sequence_list_string += sequence + ","
sequence_list_string.rstrip()

success_count = len(glob.glob(outSuccessPrefix+"*.Success"))

# test for continuity with prior run for merge runs
if continuity_test:
    if refName != old_ref_name:
        print "\nReference name is different from prior run(s)"
        print "Pipeline Stopped: please use correct reference\n"
        sys.exit()
    if refGenbank and old_ref_format != 'Genbank':
        print "\nExpected FASTA reference, not Genbank"
        print "Pipeline Stopped: please use correct reference\n"
        sys.exit()
    if not refGenbank and old_ref_format == 'Genbank':
        print "\nExpected Genbank reference, not FASTA"
        print "Pipeline Stopped: please use correct reference\n"
        sys.exit()
    if len(replicons) != int(old_replicon_count):
        print "\nDifferent number of replicons found in reference"
        print "Pipeline Stopped: please use correct reference\n"
        sys.exit()
    if runType != old_run_type:
        print "\nExpected " + old_run_type + " run; " + runType + " run specified"
        print "Pipeline Stopped: please specify same run type as previously used ("+old_run_type+")\n"
        sys.exit()
    replicon_fail = False
    if runType == 'pangenome':
        for replicon in old_core_replicon_list:
            if replicon not in core_replicons:
                replicon_fail = True
        for replicon in core_replicons:
            if replicon not in old_core_replicon_list:
                replicon_fail = True
        if replicion_fail:
            print "\nCore replicons have changed since last run"
            print "Pipeline Stopped: please specify reference with same replicons as previously used\n"
            sys.exit()
    else:
        replicon_names = []
        for replicon in replicons:
            replicon_names.append(replicon[0])
            if replicon[0] not in old_replicon_list:
                replicon_fail = True
        for replicon in old_replicon_list:
            if replicon not in replicon_names:
                replicon_fail = True
        if replicon_fail:
            print "\nReplicons have changed since last run"
            print "Pipeline Stopped: please specify reference with same replicons as previously\n"
            sys.exit()

    if old_mapped_fail != 'off':
        old_check_reads_mapped = ""
        for replicon in old_replicon_test_list:
            if len(old_replicon_test_list) > 1:
                old_check_reads_mapped += (replicon + ",")
            else:
                old_check_reads_mapped = replicon
        if len(old_replicon_test_list) > 1:
            old_check_reads_mapped += 'x,'
            for number in range(len(old_replicon_percent_list)-1):                
                if number < len(old_replicon_percent_list)-2:
                    old_check_reads_mapped += str(float(old_replicon_percent_list[number])/100) + ','
                else:
                    old_check_reads_mapped += str(float(old_replicon_percent_list[number])/100)

        if check_reads_mapped != 'off':
            check_test = True
            if len(old_replicon_test_list) == 1:
                if check_reads_mapped.find(',x,') != -1:
                    check_test = False
                else:
                    if old_replicon_test_list[0] != check_reads_mapped:
                        check_test = False

            if len(old_replicon_test_list) > 1:
                if check_reads_mapped.find(',x,') == -1:
                    check_test = False
                else:
                    all_values = check_reads_mapped.split(',x,')
                    names = all_values[0].split(',')
                    values = all_values[1].split(',')
                    if len(names) != len(old_replicon_test_list):
                        check_test = False
                    else:
                        for name in names:
                            if name not in old_replicon_test_list:
                                check_test = False
                        for name in old_replicon_test_list:
                            if name not in names:
                                check_test = False
                        replicon_value_test = []
                        for number in range(len(old_replicon_percent_list)-1):
                            replicon_value_test.append(float(old_replicon_percent_list[number])/100)
                        if len(values) != len(replicon_value_test):
                            check_test = False
                        else:
                            for number in range(len(values)):
                                if float(values[number]) != replicon_value_test[number]:
                                    check_test = False 

            if not check_test:
                print "\n'check_reads_mapped' has changed since last run"
                value_change = False
                attempts_count = 0
                while value_change == False:
                    keyboard_entry = raw_input("\nWhich 'check_reads_mapped' do you want to use: \n[1] the new value, or\n[2] the old value?\n")
                    if keyboard_entry == '1' or keyboard_entry == '2':
                        value_change = True
                    if keyboard_entry == '2':
                        print "\n'check_reads_mapped' restored to previously used setting:"
                        print old_check_reads_mapped + "\n\n"
                        check_reads_mapped = old_check_reads_mapped
                    elif keyboard_entry == '1':
                        print "\n'check_reads_mapped' using new setting:"
                        print check_reads_mapped + "\n"
                else:
                    attempts_count +=1
                    if attempts_count >= 3:
                        print "\nPipeline Stopped: too many tries\n"
                        sys.exit()
                    else:
                        print "Please enter either '1' for 'yes', or '2' for 'no'"

        else:
            print "\n'check_reads_mapped' has changed to 'off' since last run"    
            value_change = False
            attempts_count = 0
            while value_change == False:
                keyboard_entry = raw_input('\nAre you sure you want to turn off checking percentage of reads mapped: \n[1] yes, or\n[2] no?\n')
                if keyboard_entry == '1' or keyboard_entry == '2':
                    value_change = True
                    if keyboard_entry == '2':
                        print "\n'check_reads_mapped' restored to previously used setting:"
                        print old_check_reads_mapped + "\n"
                        check_reads_mapped = old_check_reads_mapped
                    elif keyboard_entry == '1':
                        print "\n'check_reads_mapped' set to 'off' confirmed\n"
                else:
                    attempts_count +=1
                    if attempts_count >= 3:
                        print "\nPipeline Stopped: too many tries\n"
                        sys.exit()
                    else:
                        print "Please enter either '1' for 'yes', or '2' for 'no'"

    else:
        if check_reads_mapped != "off":
            print "\n'check_reads_mapped' has changed from 'off' since last run"    
            value_change = False
            attempts_count = 0
            while value_change == False:
                keyboard_entry = raw_input('\nAre you sure you want to turn on checking percentage of reads mapped: \n[1] yes, or\n[2] no?\n')
                if keyboard_entry == '1' or keyboard_entry == '2':
                    value_change = True
                    if keyboard_entry == '2':
                        print "\n'check_reads_mapped' set to 'off' confirmed\n"
                        check_reads_mapped = 'off'
                    elif keyboard_entry == '1':
                        print "\n'check_reads_mapped' changed from 'off' and set to:"
                        print check_reads_mapped + "\n"
                else:
                    attempts_count +=1
                    if attempts_count >= 3:
                        print "\nPipeline Stopped: too many tries\n"
                        sys.exit()
                    else:
                        print "Please enter either '1' for 'yes', or '2' for 'no'"

    if old_bowtie_preset != '' and old_bowtie_preset != bowtie_map_type:
        print "\n'bowtie_map_type' has changed since last run"
        value_change = False
        attempts_count = 0
        while value_change == False:
            keyboard_entry = raw_input('\nDo you wish to use: \n[1] the new value ('+bowtie_map_type +'), or\n[2] the old value? ('+old_bowtie_preset+')\nNote: this only affects new read sets\n')
            if keyboard_entry == '1' or keyboard_entry == '2':
                value_change = True
                if keyboard_entry == '2':
                    bowtie_map_type = old_bowtie_preset 
                print "\n'bowtie_map_type' set to " + bowtie_map_type
            else:
                attempts_count +=1
                if attempts_count >= 3:
                    print "\nPipeline Stopped: too many tries\n"
                    sys.exit()
                else:
                    print "Please enter either '1' for new value, or '2' for old value"

    if old_bowtie_X != '' and int(old_bowtie_X) != bowtie_X_value:
        print "\n'bowtie_X_value' has changed since last run"
        value_change = False
        attempts_count = 0
        while value_change == False:
            keyboard_entry = raw_input('\nDo you wish to use: \n[1] the new value ('+str(bowtie_X_value) +'), or\n[2] the old value? ('+old_bowtie_X+')\nNote: this only affects new read sets\n')
            if keyboard_entry == '1' or keyboard_entry == '2':
                value_change = True
                if keyboard_entry == '2':
                    bowtie_X_value = int(old_bowtie_X)
                print "\n'bowtie_X_value' set to " + str(bowtie_X_value)
            else:
                attempts_count +=1
                if attempts_count >= 3:
                    print "\nPipeline Stopped: too many tries\n"
                    sys.exit()
                else:
                    print "Please enter either '1' for new value, or '2' for old value"

    if int(old_min_depth) != minDepth:
        print "\n'minimum_depth' for SNP calling has changed since last run"
        value_change = False
        attempts_count = 0
        while value_change == False:
            keyboard_entry = raw_input('\nDo you wish to use: \n[1] the new value ('+str(minDepth)+'), or\n[2] the old value? ('+old_min_depth+')\nNote: this only affects new read sets\n')
            if keyboard_entry == '1' or keyboard_entry == '2':
                value_change = True
                if keyboard_entry == '2':
                    minDepth = int(old_min_depth)
                print "\n'minimum_depth' set to " + str(minDepth)
            else:
                attempts_count +=1
                if attempts_count >= 3:
                    print "\nPipeline Stopped: too many tries\n"
                    sys.exit()
                else:
                    print "Please enter either '1' for new value, or '2' for old value"

    if int(old_cover_fail) != coverFail:
        print "\n'cover_fail' has changed since last run"
        value_change = False
        attempts_count = 0
        while value_change == False:
            keyboard_entry = raw_input('\nDo you wish to use: \n[1] the new value ('+str(coverFail)+'), or\n[2] the old value? ('+old_cover_fail+')\n')
            if keyboard_entry == '1' or keyboard_entry == '2':
                value_change = True
                if keyboard_entry == '2':
                    coverFail = int(old_cover_fail)
                print "\n'cover_fail' set to " + str(coverFail)
            else:
                attempts_count +=1
                if attempts_count >= 3:
                    print "\nPipeline Stopped: too many tries\n"
                    sys.exit()
                else:
                    print "Please enter either '1' for new value, or '2' for old value"

    if int(old_depth_fail) != depthFail:
        print "\n'depth_fail' has changed since last run"
        value_change = False
        attempts_count = 0
        while value_change == False:
            keyboard_entry = raw_input('\nDo you wish to use: \n[1] the new value ('+str(depthFail)+'), or\n[2] the old value? ('+old_depth_fail+')\n')
            if keyboard_entry == '1' or keyboard_entry == '2':
                value_change = True
                if keyboard_entry == '2':
                    depthFail = int(old_depth_fail)
                print "\n'depth_fail' set to " + str(depthFail)
            else:
                attempts_count +=1
                if attempts_count >= 3:
                    print "\nPipeline Stopped: too many tries\n"
                    sys.exit()
                else:
                    print "Please enter '1' for new value, or '2' for old value"

    if old_mapped_fail != 'off' and check_reads_mapped != 'off':
        if int(old_mapped_fail) != mappedFail:
            print "\n'mapped_fail' has changed since last run"
            value_change = False
            attempts_count = 0
            while value_change == False:
                keyboard_entry = raw_input('\nDo you wish to use: \n[1] the new value ('+str(mappedFail)+'), or\n[2] the old value? ('+old_mapped_fail+')\n')
                if keyboard_entry == '1' or keyboard_entry == '2':
                    value_change = True
                    if keyboard_entry == '2':
                        mappedFail = int(old_mapped_fail)
                    print "\n'mapped_fail' set to " + str(mappedFail)
                else:
                    attempts_count +=1
                    if attempts_count >= 3:
                        print "\nPipeline Stopped: too many tries\n"
                        sys.exit()
                    else:
                        print "Please enter either '1' for new value, or '2' for old value"

    if int(old_sd_out) != sdOutgroupMultiplier:
        print "\n'sd_out' has changed since last run"
        value_change = False
        attempts_count = 0
        while value_change == False:
            keyboard_entry = raw_input('\nDo you wish to use: \n[1] the new value ('+str(sdOutgroupMultiplier)+'), or\n[2] the old value? ('+old_sd_out+')\n')
            if keyboard_entry == '1' or keyboard_entry == '2':
                value_change = True
                if keyboard_entry == '2':
                    sdOutgroupMultiplier = old_sd_out
                print "\n'sd_out' set to " + str(sdOutgroupMultiplier)
            else:
                attempts_count +=1
                if attempts_count >= 3:
                    print "\nPipeline Stopped: too many tries\n"
                    sys.exit()
                else:
                    print "Please enter '1' for new value, or '2' for old value"

    if float(old_conservation) != conservation:
        print "\n'conservation' has changed since last run"
        value_change = False
        attempts_count = 0
        while value_change == False:
            keyboard_entry = raw_input('\nDo you wish to use: \n[1] the new value ('+str(conservation)+'), or\n[2] the old value? ('+old_conservation+')\n')
            if keyboard_entry == '1' or keyboard_entry == '2':
                value_change = True
                if keyboard_entry == '2':
                    conservation = int(old_conservation)
                print "\n'conservation' set to " + str(conservation)
            else:
                attempts_count +=1
                if attempts_count >= 3:
                    print "\nPipeline Stopped: too many tries\n"
                    sys.exit()
                else:
                    print "Please enter '1' for new value, or '2' for old value"



#Phew! Now that's all set up, we can begin...
#but first, output run conditions to user and get confirmation to run
print "\nRedDog " + version + " - " + runType + " run\n"
#print license information
print "Mapping: " + mapping_out
if mapping == 'bowtie':
    print "Preset Option: " + bowtie_map_type
if refGenbank == False:
    ref_string = 'FASTA'
else:
    ref_string = 'GenBank'
print str(len(replicons)) + " replicon(s) in " + ref_string + " reference " + refName
if runType != '':
    if runType == 'phylogeny':
        number_string = str(len(replicons))
    else:
        number_string = str(len(core_replicons))
    print number_string + " replicon(s) to be reported"
number_string = str(len(sequence_list))
if readType == 'PE':
    print number_string + " sequence pair(s) to be mapped"
else:
    print number_string + " sequence(s) to be mapped"

print "\nOutput folder:"
print outPrefix
if outMerge != '':
    print "Remember: this output folder will be deleted at the end of the run\n"
    print "Merge new sets with the following folder ('out_merge_target'):"
    print outMerge

start_run = False
start_count = 0
while start_run == False:
    keyboard_entry = raw_input('\nStart Pipeline? (y/n) ')
    if keyboard_entry == 'y' or keyboard_entry == 'Y' or keyboard_entry == 'yes' or keyboard_entry == 'Yes' or keyboard_entry == 'YES':
        start_run = True
    elif keyboard_entry == 'n' or keyboard_entry == 'N' or keyboard_entry == 'no' or keyboard_entry == 'No' or keyboard_entry == 'NO':
        print "\nPipeline Stopped: user request\n"
        sys.exit()
    else:
        start_count +=1
        if start_count >= 3:
            print "\nPipeline Stopped: too many tries\n"
            sys.exit()
        else:
            print "Please enter 'y' (yes) or 'n' (no)"

print "\nStarting pipeline..."
stage_count = 0

# Create temp and other output subfolders
@files(input==None, outSuccessPrefix + "dir.makeDir.Success")
def makeDir(input, flagFile):
    runStageCheck('makeDir', flagFile, outPrefix, full_sequence_list_string)
stage_count += 1

if refGenbank == False:
    # Copy reference to outTemp
    newReference = outTempPrefix + refName + refExt
    @follows(makeDir)
    @files(reference, [newReference, outSuccessPrefix + refName + '.copyRef.Success'])
    def copyRef(reference, outputs):
        _output, flagFile = outputs
        runStageCheck('copyRef', flagFile, reference, outTempPrefix)
    reference = newReference
    stage_count += 1

    # Index copy of reference for SNP calling by samtools
    @follows(copyRef)
    @files(reference, [reference + '.fai', outSuccessPrefix + refName + '.indexRef.Success'])
    def indexRef(reference, outputs):
        _output, flagFile = outputs
        runStageCheck('indexRef', flagFile, reference)
    stage_count += 1

    if mapping == 'bowtie':
        # Index the reference file for bowtie2
        @follows(copyRef)
        @files(reference, [outTempPrefix + refName + '.1.bt2', outSuccessPrefix + refName + '.buildBowtieIndex.Success'])
        def buildBowtieIndex(reference, outputs):
            output, flagFile = outputs
            base = output[:-6]
            runStageCheck('buildBowtieIndex', flagFile, reference, base)
        stage_count += 1
    else:
        # Index the reference file for bwa
        @follows(copyRef)
        @files(reference, [reference + '.bwt', outSuccessPrefix + refName + '.buildBWAIndex.Success'])
        def buildBWAIndex(reference, outputs):
            _output, flagFile = outputs
            runStageCheck('buildBWAIndex', flagFile, reference)
        stage_count += 1

else:
    #make a fasta version of the reference for mapping
    genbank = reference
    reference = outTempPrefix + refName + ".fasta"
    @follows(makeDir)
    @files(genbank, [reference, outSuccessPrefix + refName + '.makeRef.Success'])
    def makeRef(genbank, outputs):        
        output, flagFile = outputs
        runStageCheck('makeRef', flagFile, genbank, output)
    stage_count += 1

    # Index constructed reference for SNP calling by samtools
    @follows(makeRef)
    @files(reference, [reference + '.fai', outSuccessPrefix + refName + '.indexRef.Success'])
    def indexRef(reference, outputs):
        _output, flagFile = outputs
        runStageCheck('indexRef', flagFile, reference)
    stage_count += 1

    if mapping == 'bowtie':
        # Index the reference file for bowtie2
        @follows(makeRef)
        @files(reference, [outTempPrefix + refName + '.1.bt2', outSuccessPrefix + refName + '.buildBowtieIndex.Success'])
        def buildBowtieIndex(reference, outputs):
            output, flagFile = outputs
            base = output[:-6]
            runStageCheck('buildBowtieIndex', flagFile, reference, base)
        stage_count += 1
    else:
        # Index the reference file for bwa
        @follows(makeRef)
        @files(reference, [reference + '.bwt', outSuccessPrefix + refName + '.buildBWAIndex.Success'])
        def buildBWAIndex(reference, outputs):
            _output, flagFile = outputs
            runStageCheck('buildBWAIndex', flagFile, reference)
        stage_count += 1

if mapping == 'bowtie':
    if readType == "PE":
        #align reads with Bowtie2 to sorted bam
        input = r'(.*)\/(.+)_1\.fastq.gz'
        extraInput = [r'\1/\2_2.fastq.gz']
        @follows(buildBowtieIndex)
        @transform(sequences, regex(input), add_inputs(extraInput), [outTempPrefix + r'\2/\2.bam', outSuccessPrefix + r'\2.alignBowtiePE.Success'])
        def alignBowtiePE(inputs, outputs):
            output, flagFile = outputs
            (prefix, name, ext) = splitPath(output)
            out = prefix + '/' + name
            seq1, [seq2] = inputs
            base = outTempPrefix + refName
            runStageCheck('alignBowtiePE', flagFile, bowtie_map_type, base, seq1, seq2, bowtie_X_value, out)
        stage_count += len(sequence_list) 

        # Index sorted BAM alignments using samtools
        @transform(alignBowtiePE, regex(r"(.*)\/(.+).bam"), [r'\1/\2.bam.bai', outSuccessPrefix + r'\2.indexBam.Success'])
        def indexBam(inputs, outputs):
            output, flagFile = outputs
            bamFile, _success = inputs
            runStageCheck('indexBam', flagFile, bamFile)
        stage_count += len(sequence_list) 

        #get bam stats by replicon
        @follows(indexBam)
        @transform(alignBowtiePE, regex(r"(.*)\/(.+).bam"), [r'\1/\2_samStats.txt', outSuccessPrefix + r'\2.getSamStats.Success'])
        def getSamStats(inputs, outputs):
            output, flagFile = outputs
            bamFile, _success = inputs
            runStageCheck('getSamStats', flagFile, bamFile, output)
        stage_count += len(sequence_list) 

# checkpoint_getSamStats
        @merge(getSamStats, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "getSamStats.checkpoint.Success"])
        def checkpoint_getSamStats(inputs,outputs):
            output, flagFile = outputs
            runStageCheck('checkpoint', flagFile, outTempPrefix, 'getSamStats')
        stage_count += 1

        #filter unmapped reads
        @follows(indexBam, checkpoint_getSamStats)
        @transform(alignBowtiePE, regex(r"(.*)\/(.+).bam"), [outBamPrefix + r'\2.bam', outSuccessPrefix + r'\2.filterUnmapped.Success'])
        def filterUnmapped(inputs, outputs):
            output, flagFile = outputs
            bamFile, _success = inputs
            (prefix, name, ext) = splitPath(output)
            out = prefix + '/' + name
            runStageCheck('filterUnmapped', flagFile, bamFile, out)
        stage_count += len(sequence_list) 

    else:
        if readType=="SE":
            #align reads with Bowtie2 to sorted bam
            @follows(buildBowtieIndex)
            @transform(sequences, regex(r"(.*)\/(.+).fastq.gz"), [outTempPrefix + r'\2/\2.bam', outSuccessPrefix + r'\2.alignBowtie.Success'])
            def alignBowtie(input, outputs):
                output, flagFile = outputs
                (prefix, name, ext) = splitPath(output)
                out = prefix + '/' + name
                base = outTempPrefix + refName
                runStageCheck('alignBowtie', flagFile, bowtie_map_type, base, input, out)
            stage_count += len(sequence_list) 
     
        else: #readType == "IT"
            #align reads with Bowtie2 to sorted bam
            @follows(buildBowtieIndex)
            @transform(sequences, regex(r"(.*)\/(.+)_in.iontor.fastq.gz"), [outTempPrefix + r'\2/\2.bam', outSuccessPrefix + r'\2.alignBowtie.Success'])
            def alignBowtie(input, outputs):
                output, flagFile = outputs
                (prefix, name, ext) = splitPath(output)
                out = prefix + '/' + name
                base = outTempPrefix + refName
                runStageCheck('alignBowtie', flagFile, bowtie_map_type, base, input, out)
            stage_count += len(sequence_list) 

        # Index sorted BAM alignments using samtools
        @transform(alignBowtie, regex(r"(.*)\/(.+).bam"), [r'\1/\2.bam.bai', outSuccessPrefix + r'\2.indexBam.Success'])
        def indexBam(inputs, outputs):
            output, flagFile = outputs
            bamFile, _success = inputs
            runStageCheck('indexBam', flagFile, bamFile)
        stage_count += len(sequence_list) 

        #get bam stats by replicon
        @follows(indexBam)
        @transform(alignBowtie, regex(r"(.*)\/(.+).bam"), [r'\1/\2_samStats.txt', outSuccessPrefix + r'\2.getSamStats.Success'])
        def getSamStats(inputs, outputs):
            output, flagFile = outputs
            bamFile, _success = inputs
            runStageCheck('getSamStats', flagFile, bamFile, output)
        stage_count += len(sequence_list) 

# checkpoint_getSamStats
        @merge(getSamStats, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "getSamStats.checkpoint.Success"])
        def checkpoint_getSamStats(inputs,outputs):
            output, flagFile = outputs
            runStageCheck('checkpoint', flagFile, outTempPrefix, 'getSamStats')
        stage_count += 1

        #filter unmapped reads
        @follows(indexBam, checkpoint_getSamStats)
        @transform(alignBowtie, regex(r"(.*)\/(.+).bam"), [outBamPrefix + r'\2.bam', outSuccessPrefix + r'\2.filterUnmapped.Success'])
        def filterUnmapped(inputs, outputs):
            output, flagFile = outputs
            bamFile, _success = inputs
            (prefix, name, ext) = splitPath(output)
            out = prefix + '/' + name
            runStageCheck('filterUnmapped', flagFile, bamFile, out)
        stage_count += len(sequence_list) 

else: # mapping = 'BWA'
    if readType == "PE":
        # Align 'forward' sequences reads to the reference genome.
        @follows(buildBWAIndex)
        @transform(sequences, regex(r'(.*)\/(.+)_1\.fastq.gz'), [outTempPrefix + r'\2/\2_1.sai', outSuccessPrefix + r"\2.alignForward.Success"])
        def alignForward(sequence, outputs):
            output, flagFile = outputs
            runStageCheck('alignSequence', flagFile, reference, sequence, output)
        stage_count += len(sequences)/2 

        # Align 'reverse' sequences reads to the reference genome.
        @follows(buildBWAIndex)
        @transform(sequences, regex(r'(.*)\/(.+)_2\.fastq.gz'), [outTempPrefix + r'\2/\2_2.sai', outSuccessPrefix + r"\2.alignReverse.Success"])
        def alignReverse(sequence, outputs):
            output, flagFile = outputs
            runStageCheck('alignSequence', flagFile, reference, sequence, output)
        stage_count += len(sequences)/2

        #align reads with BWA to sorted bam
        input = r'(.*)\/(.+)_1\.sai'
        extraInput = [r'\1/\2_2.sai']
        @follows(alignReverse)
        @transform(alignForward, regex(input), add_inputs(extraInput), [r'\1/\2.bam', outSuccessPrefix + r'\2.alignBWAPE.Success'])
        def alignBWAPE(inputs, outputs):
            output, flagFile = outputs
            [align1, _success], [align2] = inputs
            (prefix, name, ext) = splitPath(align1)
            seqPrefix = ""
            for sequence in sequences:
                if sequence.find(name) != -1 and seqPrefix == "":
                    (seqPrefix, seqName, seqExt) = splitPath(sequence)
            name = seqPrefix + "/" + name[:-2]
            seq1 = name + '_1.fastq.gz'
            seq2 = name + '_2.fastq.gz'
            (prefix, name, ext) = splitPath(output)
            out = prefix + '/' + name
            runStageCheck('alignBWAPE', flagFile, reference, align1, align2, seq1, seq2, out)
        stage_count += len(sequence_list) 

        # Index sorted BAM alignments using samtools
        @transform(alignBWAPE, regex(r"(.*)\/(.+).bam"), [r'\1/\2.bam.bai', outSuccessPrefix + r'\2.indexBam.Success'])
        def indexBam(inputs, outputs):
            output, flagFile = outputs
            bamFile, _success = inputs
            runStageCheck('indexBam', flagFile, bamFile)
        stage_count += len(sequence_list) 

        #get bam stats by replicon
        @follows(indexBam)
        @transform(alignBWAPE, regex(r"(.*)\/(.+).bam"), [r'\1/\2_samStats.txt', outSuccessPrefix + r'\2.getSamStats.Success'])
        def getSamStats(inputs, outputs):
            output, flagFile = outputs
            bamFile, _success = inputs
            runStageCheck('getSamStats', flagFile, bamFile, output)
        stage_count += len(sequence_list) 

# checkpoint_getSamStats
        @merge(getSamStats, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "getSamStats.checkpoint.Success"])
        def checkpoint_getSamStats(inputs,outputs):
            output, flagFile = outputs
            runStageCheck('checkpoint', flagFile, outTempPrefix, 'getSamStats')
        stage_count += 1

        #filter unmapped reads
        @follows(indexBam, checkpoint_getSamStats)
        @transform(alignBWAPE, regex(r"(.*)\/(.+).bam"), [outBamPrefix + r'\2.bam', outSuccessPrefix + r'\2.filterUnmapped.Success'])
        def filterUnmapped(inputs, outputs):
            output, flagFile = outputs
            bamFile, _success = inputs
            (prefix, name, ext) = splitPath(output)
            out = prefix + '/' + name
            runStageCheck('filterUnmapped', flagFile, bamFile, out)
        stage_count += len(sequence_list) 

    else:
        # Align sequence reads to the reference genome.
        @follows(buildBWAIndex)
        @transform(sequences, regex(r"(.*)\/(.+).fastq.gz"), [outTempPrefix + r'\2/\2.sai', outSuccessPrefix + r"\2.alignSequence.Success"])
        def alignSequence(sequence, outputs):
            output, flagFile = outputs
            runStageCheck('alignSequence', flagFile, reference, sequence, output)
        stage_count += len(sequences)

        #align reads with BWA to sorted bam
        @follows(buildBWAIndex)
        @transform(alignSequence, regex(r"(.*)\/(.+).sai"), [r'\1/\2.bam', outSuccessPrefix + r'\2.alignBWASE.Success'])
        def alignBWASE(inputs, outputs):
            output, flagFile = outputs
            align, _success = inputs
            (prefix, name, ext) = splitPath(align)
            seqPrefix = ""
            for sequence in sequences:
                if sequence.find(name) != -1 and seqPrefix == "":
                    (seqPrefix, seqName, seqExt) = splitPath(sequence)
            seq = seqPrefix + "/" + name + '.fastq.gz'
            (prefix, name, ext) = splitPath(output)
            out = prefix + '/' + name
            runStageCheck('alignBWASE', flagFile, reference, align, seq, out)
        stage_count += len(sequence_list) 

        # Index sorted BAM alignments using samtools
        @transform(alignBWASE, regex(r"(.*)\/(.+).bam"), [r'\1/\2.bam.bai', outSuccessPrefix + r'\2.indexBam.Success'])
        def indexBam(inputs, outputs):
            output, flagFile = outputs
            bamFile, _success = inputs
            runStageCheck('indexBam', flagFile, bamFile)
        stage_count += len(sequence_list) 

        #get bam stats by replicon
        @follows(indexBam)
        @transform(alignBWASE, regex(r"(.*)\/(.+).bam"), [r'\1/\2_samStats.txt', outSuccessPrefix + r'\2.getSamStats.Success'])
        def getSamStats(inputs, outputs):
            output, flagFile = outputs
            bamFile, _success = inputs
            runStageCheck('getSamStats', flagFile, bamFile, output)
        stage_count += len(sequence_list) 

# checkpoint_getSamStats
        @merge(getSamStats, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "getSamStats.checkpoint.Success"])
        def checkpoint_callRepSNPs(inputs,outputs):
            output, flagFile = outputs
            runStageCheck('checkpoint', flagFile, outTempPrefix, 'getSamStats')
        stage_count += 1

        #filter unmapped reads
        @follows(indexBam, checkpoint_getSamStats)
        @transform(alignBWASE, regex(r"(.*)\/(.+).bam"), [outBamPrefix + r'\2.bam', outSuccessPrefix + r'\2.filterUnmapped.Success'])
        def filterUnmapped(inputs, outputs):
            output, flagFile = outputs
            bamFile, _success = inputs
            (prefix, name, ext) = splitPath(output)
            out = prefix + '/' + name
            runStageCheck('filterUnmapped', flagFile, bamFile, out)
        stage_count += len(sequence_list) 

# Index sorted BAM alignments using samtools
@transform(filterUnmapped, regex(r"(.*)\/(.+).bam"), [r'\1/\2.bam.bai', outSuccessPrefix + r'\2.indexFilteredBam.Success'])
def indexFilteredBam(inputs, outputs):
    output, flagFile = outputs
    bamFile, _success = inputs
    runStageCheck('indexBam', flagFile, bamFile)
stage_count += len(sequence_list) 

# get consensus sequence from bam
@follows(indexFilteredBam)
@follows(indexRef)
@transform(filterUnmapped, regex(r"(.*)\/(.+).bam"), [outTempPrefix + r"\2/\2_cns.fq", outSuccessPrefix + r"\2.getConsensus.Success"])
def getConsensus(inputs, outputs):
    output, flagFile = outputs
    bamFile, _success = inputs
    runStageCheck('getConsensus', flagFile, reference, bamFile, output)
stage_count += len(sequence_list) 

# Get coverage from BAM
@follows(indexFilteredBam)
@follows(indexRef)
@transform(filterUnmapped, regex(r"(.*)\/(.+).bam"), [outTempPrefix + r"\2/\2_coverage.txt", outSuccessPrefix + r"\2.getCoverage.Success"])
def getCoverage(inputs, outputs):
    output, flagFile = outputs
    bamFile, _success = inputs
    runStageCheck('getCoverage', flagFile, bamFile, output)
stage_count += len(sequence_list) 

# Get coverage by replicon
@transform(getCoverage, regex(r"(.*)\/(.+)_coverage.txt"), [r'\1/\2_rep_cover.txt', outSuccessPrefix + r'\2.getCoverByRep.Success'])
def getCoverByRep(inputs, outputs):
    output, flagFile = outputs
    coverageFile, _success = inputs
    runStageCheck('getCoverByRep', flagFile, reference, coverageFile, output)
stage_count += len(sequence_list) 

# Derive run statistics - first, the 'all statistics' data
@follows(checkpoint_getSamStats)
@transform(getCoverByRep, regex(r"(.*)\/(.+)_rep_cover.txt"), [r'\1/\2_AllStats.txt', outSuccessPrefix + r'\2.deriveAllStats.Success'])
def deriveAllStats(inputs, outputs):
    output, flagFile = outputs
    coverFile, _success = inputs
    runStageCheck('deriveAllStats', flagFile, coverFile)
stage_count += len(sequence_list)

# Collate run statistics into tab-delimited file
# and add header to the statistics file
# example rep cover file is to establish replicon order in the header (order from reference)
# also provides a user-friendly version of the all stats file with isolates as the header (for spreadsheets) 
@merge(deriveAllStats, [outPrefix + refName + "_AllStats.txt", outSuccessPrefix + refName + ".collateAllStats.Success"])
def collateAllStats(inputs, outputs):
    output, flagFile = outputs
    example_name = sequence_list[0]
    exampleRepCover = outTempPrefix + example_name + '/' + example_name + "_rep_cover.txt"
    runStageCheck('collateAllStats', flagFile, refName, exampleRepCover, outPrefix, sequence_list_string)
stage_count += 1

if runType == "pangenome":
    #create inputs for callRepSNPs
    def snpsByCoreReplicons():
        for repliconName in core_replicons:
            for seqName in sequence_list:
                sortedBam = outBamPrefix + seqName + '.bam'
                output = outTempPrefix + seqName + '/callRepSNPs/' + seqName + '_' + repliconName + '_raw.bcf'
                flagFile = outSuccessPrefix + seqName + '_' + repliconName + '.callRepSNPs.Success'
                yield([sortedBam, [output, flagFile], repliconName])

    # Call SNPs for core replicon(s)
    @follows(indexFilteredBam)
    @follows(indexRef)
    @files(snpsByCoreReplicons)
    def callRepSNPs(input, outputs, repliconName):
        output, flagFile  = outputs
        bcf_option = SNPcaller + 'v'
        runStageCheck('callRepSNPs', flagFile, reference, input, repliconName, bcf_option, output)
    stage_count += (len(sequence_list)*len(core_replicons)) 

# checkpoint_callRepSNPs
    @merge(callRepSNPs, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "callRepSNPs.checkpoint.Success"])
    def checkpoint_callRepSNPs(inputs,outputs):
        output, flagFile = outputs
        runStageCheck('checkpoint', flagFile, outTempPrefix, 'callRepSNPs')
    stage_count += 1

    #create inputs for filter variants on Q30
    def q30FilterByCoreReplicons():
        for repliconName in core_replicons:
            for seqName in sequence_list:
                rawBCF = outTempPrefix + seqName + '/callRepSNPs/' + seqName + '_' + repliconName + '_raw.bcf'
                coverFile = outTempPrefix + seqName + '/' + seqName + '_rep_cover.txt'                        
                output = outTempPrefix + seqName + '/q30VarFilter/' + seqName + '_' + repliconName + '_raw.vcf'
                flagFile = outSuccessPrefix + seqName + '_' + repliconName + '.q30VarFilter.Success'
                yield([rawBCF, [output, flagFile], coverFile, repliconName])
    
    # Filter variants on Q30
    @follows(getCoverByRep)
    @follows(checkpoint_callRepSNPs)
    @files(q30FilterByCoreReplicons)
    def q30VarFilter(rawBCF, outputs, coverFile, repliconName):
        output, flagFile = outputs
        cover = getCover(coverFile, repliconName)
        runStageCheck('q30VarFilter', flagFile, rawBCF, minDepth, cover, output)
    stage_count += (len(sequence_list)*len(core_replicons)) 

    # Filter out simple hets
    @transform(q30VarFilter, regex(r"(.*)\/(.+)_raw.vcf"), [outVcfPrefix + r"\2_q30.vcf", outSuccessPrefix + r"\2.finalFilter.Success"])
    def finalFilter(inputs, outputs):
        vcfFile,_success = inputs
        output, flagFile = outputs
        (prefix, name, ext) = splitPath(vcfFile) 
        prefix += '/hets/'       
        runStageCheck('finalFilter', flagFile, vcfFile, output, prefix, HetsVCF)
    stage_count += (len(sequence_list)*len(core_replicons)) 

    # set up for getting vcf statistics
    def vcfStatsByCoreRep():
        for repliconName in core_replicons:
            for seqName in sequence_list:
                vcfFile = outVcfPrefix + seqName + '_' + repliconName + '_q30.vcf'
                output = outTempPrefix + seqName + '/getVCFStats/' + seqName + '_' + repliconName + '_vcf.txt'
                flagFile = outSuccessPrefix + seqName + '_' + repliconName + '.getVcfStats.Success'
                yield([vcfFile, [output, flagFile]])

    # Get the vcf statistics
    @follows(finalFilter)
    @files(vcfStatsByCoreRep)
    def getVcfStats(vcfFile, outputs):
        output, flagFile = outputs
        runStageCheck('getVcfStats', flagFile, vcfFile, output)
    stage_count += (len(sequence_list)*len(core_replicons))

    # Derive run statistics - second, the statistics for each replicon (largest or user-defined list)
    def statsByCoreRep():
        for repliconName in core_replicons:
            for seqName in sequence_list:
                coverFile = outTempPrefix + seqName + '/' + seqName + '_rep_cover.txt'
                output = outTempPrefix + seqName + '/deriveRepStats/' + seqName + '_' + repliconName + '_RepStats.txt'
                flagFile = outSuccessPrefix + seqName + '_' + repliconName + '.deriveRepStats.Success'
                yield([coverFile, [output, flagFile], repliconName, depthFail, coverFail])

    @follows(getVcfStats)
    @follows(checkpoint_getSamStats)
    @follows(getCoverByRep)
    @files(statsByCoreRep)
    def deriveRepStats(coverFile, outputs, repliconName, depthFail, coverFail):
        output, flagFile = outputs
        runStageCheck('deriveRepStats', flagFile, coverFile, repliconName, depthFail, coverFail, runType, mappedFail, check_reads_mapped)
    stage_count += (len(sequence_list)*len(core_replicons)) 

# checkpoint_deriveRepStats
    @merge(deriveRepStats, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "deriveRepStats.checkpoint.Success"])
    def checkpoint_deriveRepStats(inputs,outputs):
        output, flagFile = outputs
        runStageCheck('checkpoint', flagFile, outTempPrefix, 'deriveRepStats')
    stage_count += 1

    def inputByCoreRep():
        for repliconName in core_replicons:
            output = outPrefix + refName + '_' + repliconName + '_RepStats.txt'
            flagFile = outSuccessPrefix + refName + '_' + repliconName + '.collateRepStats.Success'
            yield([input==None, [output, flagFile], refName, repliconName])

    @follows(checkpoint_deriveRepStats)
    @files(inputByCoreRep)
    def collateRepStats(input, outputs, refName, repliconName):
        output, flagFile = outputs
        runStageCheck('collateRepStats', flagFile, refName, outPrefix, repliconName, sdOutgroupMultiplier, runType, sequence_list_string)
    stage_count += len(core_replicons) 

else: # runType == "phylogeny"
    #create inputs for callRepSNPs
    def snpsByReplicons():
        for repliconName in replicons:
            for seqName in sequence_list:
                sortedBam = outBamPrefix + seqName + '.bam'
                output = outTempPrefix + seqName + '/callRepSNPs/' + seqName + '_' + repliconName[0] + '_raw.bcf'
                flagFile = outSuccessPrefix + seqName + '_' + repliconName[0] + '.callRepSNPs.Success'
                replicon = repliconName[0]
                yield([sortedBam, [output, flagFile], replicon])

    # Call SNPs for all replicon(s)
    @follows(indexFilteredBam)
    @follows(indexRef)
    @files(snpsByReplicons)
    def callRepSNPs(input, outputs, replicon):
        output, flagFile  = outputs
        bcf_option = SNPcaller + 'v'
        runStageCheck('callRepSNPs', flagFile, reference, input, replicon, bcf_option, output)
    stage_count += (len(sequence_list)*len(replicons)) 

# checkpoint_callRepSNPs
    @merge(callRepSNPs, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "callRepSNPs.checkpoint.Success"])
    def checkpoint_callRepSNPs(inputs,outputs):
        output, flagFile = outputs
        runStageCheck('checkpoint', flagFile, outTempPrefix, 'callRepSNPs')
    stage_count += 1

    #create inputs for filter variants on Q30
    def q30FilterByReplicons():
        for repliconName in replicons:
            for seqName in sequence_list:
                rawBCF = outTempPrefix + seqName + '/callRepSNPs/' + seqName + '_' + repliconName[0] + '_raw.bcf'
                coverFile = outTempPrefix + seqName + '/' + seqName + '_rep_cover.txt'    
                output = outTempPrefix + seqName + '/q30VarFilter/' + seqName + '_' + repliconName[0] + '_raw.vcf'
                flagFile = outSuccessPrefix + seqName + '_' + repliconName[0] + '.q30VarFilter.Success'
                replicon = repliconName[0]
                yield([rawBCF, [output, flagFile], coverFile, replicon])
    
    # Filter variants on Q30
    @follows(getCoverByRep)
    @follows(checkpoint_callRepSNPs)
    @files(q30FilterByReplicons)
    def q30VarFilter(rawBCF, outputs, coverFile, replicon):
        output, flagFile = outputs
        cover = getCover(coverFile, replicon)
        runStageCheck('q30VarFilter', flagFile, rawBCF, minDepth, cover, output)
    stage_count += (len(sequence_list)*len(replicons)) 

    # Filter out simple hets
    @transform(q30VarFilter, regex(r"(.*)\/(.+)_raw.vcf"), [outVcfPrefix + r"\2_q30.vcf", outSuccessPrefix + r"\2.finalFilter.Success"])
    def finalFilter(inputs, outputs):
        output, flagFile = outputs
        vcfFile, _Success = inputs
        (prefix, name, ext) = splitPath(vcfFile) 
        prefix += '/hets/'       
        runStageCheck('finalFilter', flagFile, vcfFile, output, prefix, HetsVCF)
    stage_count += (len(sequence_list)*len(replicons)) 

    # set up for getting vcf statistics
    def vcfStatsByRep():
        for repliconName in replicons:
            for seqName in sequence_list:
                vcfFile = outVcfPrefix + seqName + '_' + repliconName[0] + '_q30.vcf'
                output = outTempPrefix + seqName + '/getVCFStats/' + seqName + '_' + repliconName[0] + '_vcf.txt'
                flagFile = outSuccessPrefix + seqName + '_' + repliconName[0] + '.getVcfStats.Success'
                yield([vcfFile, [output, flagFile]])

    # Get the vcf statistics
    @follows(finalFilter)
    @files(vcfStatsByRep)
    def getVcfStats(vcfFile, outputs):
        output, flagFile = outputs
        runStageCheck('getVcfStats', flagFile, vcfFile, output)
    stage_count += (len(sequence_list)*len(replicons))

    # Derive run statistics - second, the statistics for each replicon (largest or user-defined list)
    def statsByRep():
        for repliconName in replicons:
            for seqName in sequence_list:
                coverFile = outTempPrefix + seqName + '/' + seqName + '_rep_cover.txt'
                output = outTempPrefix + seqName + '/deriveRepStats/' + seqName + '_' + repliconName[0] + '_RepStats.txt'
                flagFile = outSuccessPrefix + seqName + '_' + repliconName[0] + '.deriveRepStats.Success'
                replicon = repliconName[0]
                yield([coverFile, [output, flagFile], replicon, depthFail, coverFail])

    @follows(getVcfStats)
    @follows(checkpoint_getSamStats)
    @follows(getCoverByRep)
    @files(statsByRep)
    def deriveRepStats(coverFile, outputs, replicon, depthFail, coverFail):
        output, flagFile = outputs
        runStageCheck('deriveRepStats', flagFile, coverFile, replicon, depthFail, coverFail, runType, mappedFail, check_reads_mapped)
    stage_count += (len(sequence_list)*len(replicons)) 

# checkpoint_deriveRepStats
    @merge(deriveRepStats, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "deriveRepStats.checkpoint.Success"])
    def checkpoint_deriveRepStats(inputs,outputs):
        output, flagFile = outputs
        runStageCheck('checkpoint', flagFile, outTempPrefix, 'deriveRepStats')
    stage_count += 1

    def inputByRep():
        for repliconName in replicons:
            output = outPrefix + refName + '_' + repliconName[0] + '_RepStats.txt'
            flagFile = outSuccessPrefix + refName + '_' + repliconName[0] + '.collateRepStats.Success'
            replicon = repliconName[0]
            yield([input==None, [output, flagFile], refName, replicon])

    @follows(checkpoint_deriveRepStats)
    @files(inputByRep)
    def collateRepStats(input, outputs, refName, replicon):
        output, flagFile = outputs
        runStageCheck('collateRepStats', flagFile, refName, outPrefix, replicon, sdOutgroupMultiplier, runType, sequence_list_string)
    stage_count += len(replicons)

if outMerge != "":
    # Merge bams and vcfs to correct output merge directory
    @follows(collateAllStats)
    @follows(collateRepStats)
    @files(input==None, outSuccessPrefix + 'mergeOutputs.Success')
    def mergeOutputs(input, flagFile):
        inputBam = outBamPrefix + '*.bam'
        inputIndex = outBamPrefix + '*.bai'
        inputVcf = outVcfPrefix + '*.vcf' 
        runStageCheck('mergeOutputs', flagFile, inputBam, outMergeBam, inputIndex, outMergeBam, inputVcf, outMergeVcf)
    stage_count += 1
 
    # Merge AllStats.txt files <- need to add 'user' version when merging
    @follows(mergeOutputs)
    @transform(collateAllStats, regex(r"(.*)\/(.+)_AllStats.txt"), [outMerge + r"\2_AllStats.txt", outSuccessPrefix + r"\2.mergeAllStats.Success"])
    def mergeAllStats(inputs, outputs):
        input, _success = inputs
        output, flagFile = outputs
        runStageCheck('mergeAllStats', flagFile, input, outMerge)
    stage_count += 1

    # Merge RepStats.txt files
    @follows(mergeOutputs)
    @transform(collateRepStats, regex(r"(.*)\/(.+)_RepStats.txt"), [outMerge + r"\2_RepStats.txt", outSuccessPrefix + r"\2.mergeRepStats.Success"])
    def mergeRepStats(inputs, outputs):
        input, _success = inputs
        output, flagFile = outputs
        runStageCheck('mergeRepStats', flagFile, input, sdOutgroupMultiplier, replaceReads, outMerge, runType)
    if runType == "phylogeny":
        stage_count += len(replicons) 
    else:
        stage_count += len(core_replicons) 

    if runType == "phylogeny":
        # Start of phylogeny analysis
        def snpListByRep():
            for repliconName in replicons:
                input  = outMerge + refName + '_' + repliconName[0] + '_RepStats.txt'
                output = outTempPrefix + 'getRepSNPList/' + refName + '_' + repliconName[0] + '_SNPList.txt'
                flagFile = outSuccessPrefix + refName + '_' + repliconName[0] + '.getRepSNPList.Success'
                replicon = repliconName[0]
                yield([input, [output, flagFile], replicon])

        # Get unique SNP list for new set using stats file for each replicon
        @follows(mergeRepStats)
        @follows(mergeAllStats)
        @files(snpListByRep)
        def getRepSNPList(input, outputs, replicon):
            output, flagFile = outputs
            runStageCheck('getRepSNPList', flagFile, input, replicon, output)
        stage_count += len(replicons) 

    else: #runType == "pangenome":
        # Start of pangenome analysis
        def snpListByCoreRep():
            for repliconName in core_replicons:
                input  = outMerge + refName + '_' + repliconName + '_RepStats.txt'
                output = outTempPrefix + 'getRepSNPList/' + refName + '_' + repliconName + '_SNPList.txt'
                flagFile = outSuccessPrefix + refName + '_' + repliconName + '.getRepSNPList.Success'
                replicon = repliconName
                yield([input, [output, flagFile], replicon])

        # Get unique SNP list for new set using stats file for each replicon
        @follows(mergeRepStats)
        @follows(mergeAllStats)
        @files(snpListByCoreRep)
        def getRepSNPList(input, outputs, replicon):
            output, flagFile = outputs
            runStageCheck('getRepSNPList', flagFile, input, replicon, output)
        stage_count += len(core_replicons) 

else: #    if mergeReads == "":
    if runType == "phylogeny":
        # Start of new run phylogeny analysis
        def snpListByRep():
            for repliconName in replicons:
                input  = outPrefix + refName + '_' + repliconName[0] + '_RepStats.txt'
                output = outTempPrefix + 'getRepSNPList/' + refName + '_' + repliconName[0] + '_SNPList.txt'
                flagFile = outSuccessPrefix + refName + '_' + repliconName[0] + '.getRepSNPList.Success'
                replicon = repliconName[0]
                yield([input, [output, flagFile], replicon])

        # Get unique SNP list for new set using stats file for each replicon
        @follows(collateRepStats)
        @follows(collateAllStats)
        @files(snpListByRep)
        def getRepSNPList(input, outputs, replicon):
            output, flagFile = outputs
            runStageCheck('getRepSNPList', flagFile, input, replicon, output)
        stage_count += len(replicons) 

    else: #runType == "pangenome":
        # Start of new run pangenome analysis
        def snpListByCoreRep():
            for repliconName in core_replicons:
                input  = outPrefix + refName + '_' + repliconName + '_RepStats.txt'
                output = outTempPrefix + 'getRepSNPList/' + refName + '_' + repliconName + '_SNPList.txt'
                flagFile = outSuccessPrefix + refName + '_' + repliconName + '.getRepSNPList.Success'
                replicon = repliconName
                yield([input, [output, flagFile], replicon])

        # Get unique SNP list for new set using stats file for each replicon
        @follows(collateRepStats)
        @follows(collateAllStats)
        @files(snpListByCoreRep)
        def getRepSNPList(input, outputs, replicon):
            output, flagFile = outputs
            runStageCheck('getRepSNPList', flagFile, input, replicon, output)
        stage_count += len(core_replicons) 

# checkpoint_getRepSNPList
@merge(getRepSNPList, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "getRepSNPList.checkpoint.Success"])
def checkpoint_getRepSNPList(inputs,outputs):
    output, flagFile = outputs
    runStageCheck('checkpoint', flagFile, outTempPrefix, 'getRepSNPList')
stage_count += 1

if refGenbank == True:
    if outMerge != "":

        bams = []        
        for sequence in full_sequence_list:
            if sequence not in sequence_list:
                bams.append(outMergeBam + sequence + '.bam')

        # generate data for the gene cover and depth matrices
        @transform(getCoverage, regex(r"(.*)\/(.+)_coverage.txt"), [outTempPrefix + r'\2/\2_CoverDepthMatrix.txt', outSuccessPrefix + r'\2.deriveAllRepGeneCover.Success'])
        def deriveAllRepGeneCover(inputs, outputs):
            output, flagFile = outputs
            input, _success = inputs
            runStageCheck('deriveAllRepGeneCover', flagFile, outTempPrefix, genbank, input)
        stage_count += len(sequence_list)

        # merge the data with the 'old' gene cover and depth matrices
        @merge(deriveAllRepGeneCover, [outMerge + refName + "_CoverMatrix.csv", outSuccessPrefix + refName + "_CoverMatrix.mergeAllRepGeneCover.Success"])
        def mergeAllRepGeneCover(inputs, outputs):
            output, flagFile = outputs
            runStageCheck('mergeAllRepGeneCover', flagFile, outTempPrefix, outMerge, refName, sequence_list_string)
        stage_count += 1

        # parse gene cover (and depth - to come) matrices to summarise gene content
        @transform(mergeAllRepGeneCover, regex(r"(.*)\/(.+)_CoverMatrix.csv"), [outMerge + r"\2_GeneSummary.csv", outSuccessPrefix + r"\2_alleles.parseGeneContent.Success"])
        def parseGeneContent(inputs, outputs):
            output, flagFile = outputs
            (prefix, name, ext) = splitPath(output)
            output2 = prefix + "/" + name[:-11] + "PresenceAbsence.csv"
            input, _success = inputs
            runStageCheck('parseGeneContent', flagFile, input, output, output2)
        stage_count += 1

        # get consensus sequences for merged set
        @follows(checkpoint_getRepSNPList)
        @transform(bams, regex(r"(.*)\/(.+).bam"), [outTempPrefix + r"\2/\2_cns.fq", outSuccessPrefix + r"\2.getMergeConsensus.Success"])
        def getMergeConsensus(input, outputs):
            output, flagFile = outputs
            runStageCheck('getConsensus', flagFile, reference, input, output)
        stage_count += len(bams)

# checkpoint_getMergeConsensus
        @merge(getMergeConsensus, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "getMergeConsensus.checkpoint.Success"])
        def checkpoint_getMergeConsensus(inputs,outputs):
            output, flagFile = outputs
            runStageCheck('checkpoint', flagFile, outTempPrefix, 'getMergeConsensus')
        stage_count += 1
 
        if runType == 'pangenome':
            # generate the allele matrix entry for each isolate for each replicon
            def matrixEntryByCoreRep():
                for isolate in full_sequence_list:
                    for repliconName in core_replicons:
                        input = outTempPrefix + 'getRepSNPList/' + refName + '_' + repliconName + '_SNPList.txt'
                        output  = outTempPrefix + isolate + '/deriveRepAlleleMartix/' + refName + '_' + repliconName + '_' + isolate + '_alleles.txt'
                        flagFile = outSuccessPrefix + refName + '_' + repliconName + '_' + isolate + '.deriveRepAlleleMatrix.Success'
                        replicon = repliconName
                        consensus =  outTempPrefix + isolate + '/' + isolate + '_cns.fq'
                        repliconStats = outMerge + refName + '_' + repliconName + '_RepStats.txt'
                        merge_prefix = outMerge
                        yield([input, [output, flagFile], replicon, consensus,repliconStats, merge_prefix])

            @follows(getConsensus, checkpoint_getMergeConsensus)
            @follows(checkpoint_getRepSNPList)
            @files(matrixEntryByCoreRep)
            def deriveRepAlleleMatrix(input, outputs, replicon, consensus, repliconStats, merge_prefix):
                output, flagFile = outputs
                runStageCheck('deriveRepAlleleMatrix', flagFile, input, output, reference, replicon, consensus, repliconStats, merge_prefix)
            stage_count += (len(full_sequence_list)*len(core_replicons))

# checkpoint_deriveRepAlleleMatrix
            @merge(deriveRepAlleleMatrix, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "deriveRepAlleleMatrix.checkpoint.Success"])
            def checkpoint_deriveRepAlleleMatrix(inputs,outputs):
                output, flagFile = outputs
                runStageCheck('checkpoint', flagFile, outTempPrefix, 'deriveRepAlleleMatrix')
            stage_count += 1

            def matrixByCoreRep():
                for repliconName in core_replicons:
                    output  = outTempPrefix + refName + '_' + repliconName + '_alleles.csv'
                    flagFile = outSuccessPrefix + refName + '_' + repliconName + '.collateRepAlleleMatrix.Success'
                    rep_name = refName + '_' + repliconName
                    yield([input==None, [output, flagFile], rep_name])

            @follows(checkpoint_deriveRepAlleleMatrix)
            @files(matrixByCoreRep)
            def collateRepAlleleMatrix(input, outputs, rep_name):
                output, flagFile = outputs
                runStageCheck('collateRepAlleleMatrix', flagFile, outTempPrefix, output, full_sequence_list_string, rep_name)
            stage_count += len(core_replicons)

        else: #runType == 'phylogeny'
            # generate the allele matrix entry for each isolate for each replicon
            def matrixEntryByRep():
                for isolate in full_sequence_list:
                    for repliconName in replicons:
                        input = outTempPrefix + 'getRepSNPList/' + refName + '_' + repliconName[0] + '_SNPList.txt'
                        output  = outTempPrefix + isolate + '/deriveRepAlleleMartix/' + refName + '_' + repliconName[0] + '_' + isolate + '_alleles.txt'
                        flagFile = outSuccessPrefix + refName + '_' + repliconName[0] + '_' + isolate + '.deriveRepAlleleMatrix.Success'
                        replicon = repliconName[0]
                        consensus =  outTempPrefix + isolate + '/' + isolate + '_cns.fq'
                        repliconStats = outMerge + refName + '_' + repliconName[0] + '_RepStats.txt'
                        merge_prefix = outMerge
                        yield([input, [output, flagFile], replicon, consensus,repliconStats, merge_prefix])

            @follows(getConsensus, checkpoint_getMergeConsensus)
            @follows(checkpoint_getRepSNPList)
            @files(matrixEntryByRep)
            def deriveRepAlleleMatrix(input, outputs, replicon, consensus, repliconStats, merge_prefix):
                output, flagFile = outputs
                runStageCheck('deriveRepAlleleMatrix', flagFile, input, output, reference, replicon, consensus, repliconStats, merge_prefix)
            stage_count += (len(full_sequence_list)*len(replicons))

# checkpoint_deriveRepAlleleMatrix
            @merge(deriveRepAlleleMatrix, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "deriveRepAlleleMatrix.checkpoint.Success"])
            def checkpoint_deriveRepAlleleMatrix(inputs,outputs):
                output, flagFile = outputs
                runStageCheck('checkpoint', flagFile, outTempPrefix, 'deriveRepAlleleMatrix')
            stage_count += 1

            def matrixByRep():
                for repliconName in replicons:
                    output  = outTempPrefix + refName + '_' + repliconName[0] + '_alleles.csv'
                    flagFile = outSuccessPrefix + refName + '_' + repliconName[0] + '.collateRepAlleleMatrix.Success'
                    rep_name = refName + '_' + repliconName[0]
                    yield([input==None, [output, flagFile], rep_name])

            @follows(checkpoint_deriveRepAlleleMatrix)
            @files(matrixByRep)
            def collateRepAlleleMatrix(input, outputs, rep_name):
                output, flagFile = outputs
                runStageCheck('collateRepAlleleMatrix', flagFile, outTempPrefix, output, full_sequence_list_string, rep_name)
            stage_count += len(replicons)

        # parse SNP table to create alignment for tree and get coding consequences for snps
        @transform(collateRepAlleleMatrix, regex(r"(.*)\/(.+)_alleles.csv"), [outMerge + r"\2_alleles_var_cons"+str(conservation)+".csv", outSuccessPrefix + r"\2_alleles.parseSNPs.Success"])
        def parseSNPs(inputs, outputs):
            output, flagFile = outputs
            input,_success = inputs
            (prefix, name, ext) = splitPath(input)
            replicon = name[len(refName)+1:-8]
            runStageCheck('parseSNPs', flagFile, input, str(conservation), genbank, replicon, outMerge)
        if runType == "phylogeny":
            stage_count += len(replicons)
        else:
            stage_count += len(core_replicons)

        if conservation != 0.95:
            @transform(collateRepAlleleMatrix, regex(r"(.*)\/(.+)_alleles.csv"), [outMerge + r"\2_alleles_var_cons" + "0.95.csv", outSuccessPrefix + r"\2_alleles.parseSNPs_95.Success"])
            def parseSNPs_95(inputs, outputs):
                output, flagFile = outputs
                input,_success = inputs
                (prefix, name, ext) = splitPath(input)
                replicon = name[len(refName)+1:-8]
                conservation_temp = 0.95
                runStageCheck('parseSNPs', flagFile, input, str(conservation_temp), genbank, replicon, outMerge)
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)
            
            # create distance matrices based on pair-wise differences in SNPs
            if DifferenceMatrix:
                @follows(parseSNPs_95)
                @transform(parseSNPs, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outMerge + r"\2_SNP_diff.nxs", outSuccessPrefix + r"\2_alleles.getDifferenceMatrix.Success"])        
                def getDifferenceMatrix(inputs, outputs):
                    output, flagFile = outputs
                    input, _success = inputs
                    runStageCheck('getDifferenceMatrix', flagFile, input)
                if runType == "phylogeny":
                    stage_count += len(replicons)
                else:
                    stage_count += len(core_replicons)

        elif DifferenceMatrix:        
            # create distance matrices based on pair-wise differences in SNPs
            @transform(parseSNPs, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outMerge + r"\2_SNP_diff.nxs", outSuccessPrefix + r"\2_alleles.getDifferenceMatrix.Success"])        
            def getDifferenceMatrix(inputs, outputs):
                output, flagFile = outputs
                input, _success = inputs
                runStageCheck('getDifferenceMatrix', flagFile, input)        
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)

        # generate tree
        if conservation != 0.95:
            @follows(parseSNPs_95)
            @transform(parseSNPs, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outMerge + r"\2_alleles_var_cons"+str(conservation)+".tree", outSuccessPrefix + r"\2_alleles.makeTree.Success"])
            def makeTree(inputs, outputs):
                output, flagFile = outputs
                input, _success = inputs
                input = input[:-4] + ".mfasta"
                runStageCheck('makeTree', flagFile, input, output)
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)
        else:
            @transform(parseSNPs, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outMerge + r"\2_alleles_var_cons"+str(conservation)+".tree", outSuccessPrefix + r"\2_alleles.makeTree.Success"])
            def makeTree(inputs, outputs):
                output, flagFile = outputs
                input, _success = inputs
                input = input[:-4] + ".mfasta"
                runStageCheck('makeTree', flagFile, input, output)
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)

    else: #ie. mergeReads == "" and refGenbank == True
        # generate the gene cover and depth matrices
        @transform(getCoverage, regex(r"(.*)\/(.+)_coverage.txt"), [outTempPrefix + r'\2/\2_CoverDepthMatrix.txt', outSuccessPrefix + r'\2.deriveAllRepGeneCover.Success'])
        def deriveAllRepGeneCover(inputs, outputs):
            output, flagFile = outputs
            input, _success = inputs
            runStageCheck('deriveAllRepGeneCover', flagFile, outTempPrefix, genbank, input)
        stage_count += len(sequence_list)

        @merge(deriveAllRepGeneCover, [outPrefix + refName + "_CoverMatrix.csv", outSuccessPrefix + refName + "_CoverMatrix.collateAllRepGeneCover.Success"])
        def collateAllRepGeneCover(inputs, outputs):
            output, flagFile = outputs
            runStageCheck('collateAllRepGeneCover', flagFile, outTempPrefix, outPrefix, refName, sequence_list_string)
        stage_count += 1

        # parse gene cover (and depth - to come) matrices to summarise gene content
        @transform(collateAllRepGeneCover, regex(r"(.*)\/(.+)_CoverMatrix.csv"), [outPrefix + r"\2_GeneSummary.csv", outSuccessPrefix + r"\2_alleles.parseGeneContent.Success"])
        def parseGeneContent(inputs, outputs):
            output, flagFile = outputs
            (prefix, name, ext) = splitPath(output)
            output2 = prefix + "/" + name[:-11] + "PresenceAbsence.csv"
            input, _success = inputs
            runStageCheck('parseGeneContent', flagFile, input, output, output2)
        stage_count += 1

        if runType == 'pangenome':
            # generate the allele matrix entry for each isolate for each replicon
            def matrixEntryByCoreRep():
                for isolate in full_sequence_list:
                    for repliconName in core_replicons:
                        input = outTempPrefix + 'getRepSNPList/' + refName + '_' + repliconName + '_SNPList.txt'
                        output  = outTempPrefix + isolate + '/deriveRepAlleleMartix/' + refName + '_' + repliconName + '_' + isolate + '_alleles.txt'
                        flagFile = outSuccessPrefix + refName + '_' + repliconName + '_' + isolate + '.deriveRepAlleleMatrix.Success'
                        replicon = repliconName
                        consensus =  outTempPrefix + isolate + '/' + isolate + '_cns.fq'
                        repliconStats = outPrefix + refName + '_' + repliconName + '_RepStats.txt'
                        merge_prefix = '-'
                        yield([input, [output, flagFile], replicon, consensus,repliconStats, merge_prefix])

            @follows(getConsensus)
            @follows(checkpoint_getRepSNPList)
            @files(matrixEntryByCoreRep)
            def deriveRepAlleleMatrix(input, outputs, replicon, consensus, repliconStats, merge_prefix):
                output, flagFile = outputs
                runStageCheck('deriveRepAlleleMatrix', flagFile, input, output, reference, replicon, consensus, repliconStats, merge_prefix)
            stage_count += (len(full_sequence_list)*len(core_replicons))

# checkpoint_deriveRepAlleleMatrix
            @merge(deriveRepAlleleMatrix, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "deriveRepAlleleMatrix.checkpoint.Success"])
            def checkpoint_deriveRepAlleleMatrix(inputs,outputs):
                output, flagFile = outputs
                runStageCheck('checkpoint', flagFile, outTempPrefix, 'deriveRepAlleleMatrix')
            stage_count += 1

            def matrixByCoreRep():
                for repliconName in core_replicons:
                    output  = outTempPrefix + refName + '_' + repliconName + '_alleles.csv'
                    flagFile = outSuccessPrefix + refName + '_' + repliconName + '.collateRepAlleleMatrix.Success'
                    rep_name = refName + '_' + repliconName
                    yield([input==None, [output, flagFile], rep_name])

            @follows(checkpoint_deriveRepAlleleMatrix)
            @files(matrixByCoreRep)
            def collateRepAlleleMatrix(input, outputs, rep_name):
                output, flagFile = outputs
                runStageCheck('collateRepAlleleMatrix', flagFile, outTempPrefix, output, full_sequence_list_string, rep_name)
            stage_count += len(core_replicons)

        else: #runType == 'phylogeny'
            # generate the allele matrix entry for each isolate for each replicon
            def matrixEntryByRep():
                for isolate in full_sequence_list:
                    for repliconName in replicons:
                        input = outTempPrefix + 'getRepSNPList/' + refName + '_' + repliconName[0] + '_SNPList.txt'
                        output  = outTempPrefix + isolate + '/deriveRepAlleleMartix/' + refName + '_' + repliconName[0] + '_' + isolate + '_alleles.txt'
                        flagFile = outSuccessPrefix + refName + '_' + repliconName[0] + '_' + isolate + '.deriveRepAlleleMatrix.Success'
                        replicon = repliconName[0]
                        consensus =  outTempPrefix + isolate + '/' + isolate + '_cns.fq'
                        repliconStats = outPrefix + refName + '_' + repliconName[0] + '_RepStats.txt'
                        merge_prefix = '-'
                        yield([input, [output, flagFile], replicon, consensus,repliconStats, merge_prefix])

            @follows(getConsensus)
            @follows(checkpoint_getRepSNPList)
            @files(matrixEntryByRep)
            def deriveRepAlleleMatrix(input, outputs, replicon, consensus, repliconStats, merge_prefix):
                output, flagFile = outputs
                runStageCheck('deriveRepAlleleMatrix', flagFile, input, output, reference, replicon, consensus, repliconStats, merge_prefix)
            stage_count += (len(full_sequence_list)*len(replicons))

# checkpoint_deriveRepAlleleMatrix
            @merge(deriveRepAlleleMatrix, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "deriveRepAlleleMatrix.checkpoint.Success"])
            def checkpoint_deriveRepAlleleMatrix(inputs,outputs):
                output, flagFile = outputs
                runStageCheck('checkpoint', flagFile, outTempPrefix, 'deriveRepAlleleMatrix')
            stage_count += 1

            def matrixByRep():
                for repliconName in replicons:
                    output  = outTempPrefix + refName + '_' + repliconName[0] + '_alleles.csv'
                    flagFile = outSuccessPrefix + refName + '_' + repliconName[0] + '.collateRepAlleleMatrix.Success'
                    rep_name = refName + '_' + repliconName[0]
                    yield([input==None, [output, flagFile], rep_name])

            @follows(checkpoint_deriveRepAlleleMatrix)
            @files(matrixByRep)
            def collateRepAlleleMatrix(input, outputs, rep_name):
                output, flagFile = outputs
                runStageCheck('collateRepAlleleMatrix', flagFile, outTempPrefix, output, full_sequence_list_string, rep_name)
            stage_count += len(replicons)

        # parse SNP table to create alignment for tree and get coding consequences for snps
        @transform(collateRepAlleleMatrix, regex(r"(.*)\/(.+)_alleles.csv"), [outPrefix + r"\2_alleles_var_cons"+str(conservation)+".csv", outSuccessPrefix + r"\2_alleles.parseSNPs.Success"])
        def parseSNPs(inputs, outputs):
            output, flagFile = outputs
            input,_success = inputs
            (prefix, name, ext) = splitPath(input)
            replicon = name[len(refName)+1:-8]
            runStageCheck('parseSNPs', flagFile, input, str(conservation), genbank, replicon, outPrefix)
        if runType == "phylogeny":
            stage_count += len(replicons)
        else:
            stage_count += len(core_replicons)

        if conservation != 0.95:
            @transform(collateRepAlleleMatrix, regex(r"(.*)\/(.+)_alleles.csv"), [outPrefix + r"\2_alleles_var_cons0.95.csv", outSuccessPrefix + r"\2_alleles.parseSNPs_95.Success"])
            def parseSNPs_95(inputs, outputs):
                output, flagFile = outputs
                input,_success = inputs
                (prefix, name, ext) = splitPath(input)
                replicon = name[len(refName)+1:-8]
                conservation_temp = 1.0
                runStageCheck('parseSNPs', flagFile, input, str(conservation_temp), genbank, replicon, outPrefix)
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)
            
            # create distance matrices based on pair-wise differences in SNPs
            if DifferenceMatrix:
                @follows(parseSNPs_95)
                @transform(parseSNPs, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outPrefix + r"\2_SNP_diff.nxs", outSuccessPrefix + r"\2_alleles.getDifferenceMatrix.Success"])        
                def getDifferenceMatrix(inputs, outputs):
                    output, flagFile = outputs
                    input, _success = inputs
                    runStageCheck('getDifferenceMatrix', flagFile, input)
                if runType == "phylogeny":
                    stage_count += len(replicons)
                else:
                    stage_count += len(core_replicons)

        elif DifferenceMatrix:        
            # create distance matrices based on pair-wise differences in SNPs
            @transform(parseSNPs, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outPrefix + r"\2_SNP_diff.nxs", outSuccessPrefix + r"\2_alleles.getDifferenceMatrix.Success"])        
            def getDifferenceMatrix(inputs, outputs):
                output, flagFile = outputs
                input, _success = inputs
                runStageCheck('getDifferenceMatrix', flagFile, input)        
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)

        if conservation != 0.95:
        # generate tree
            @follows(parseSNPs_95)
            @transform(parseSNPs, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outPrefix + r"\2_alleles_var_cons"+str(conservation)+".tree", outSuccessPrefix + r"\2_alleles.makeTree.Success"])
            def makeTree(inputs, outputs):
                output, flagFile = outputs
                input, _success = inputs
                input = input[:-4] + ".mfasta"
                runStageCheck('makeTree', flagFile, input, output)
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)
        else:
            @transform(parseSNPs, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outPrefix + r"\2_alleles_var_cons"+str(conservation)+".tree", outSuccessPrefix + r"\2_alleles.makeTree.Success"])
            def makeTree(inputs, outputs):
                output, flagFile = outputs
                input, _success = inputs
                input = input[:-4] + ".mfasta"
                runStageCheck('makeTree', flagFile, input, output)
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)

else: # refGenbank == False
    if outMerge != "":
        # get consensus sequence for merged set
        bams = []
        for sequence in full_sequence_list:
            if sequence not in sequence_list:
                bams.append(outMergeBam + sequence + '.bam')

        @follows(checkpoint_getRepSNPList)
        @transform(bams, regex(r"(.*)\/(.+).bam"), [outTempPrefix + r"\2/\2_cns.fq", outSuccessPrefix + r"\2.getMergeConsensus.Success"])
        def getMergeConsensus(input, outputs):
            output, flagFile = outputs
            runStageCheck('getConsensus', flagFile, reference, input, output)
        stage_count += len(bams)

# checkpoint_getMergeConsensus
        @merge(getMergeConsensus, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "getMergeConsensus.checkpoint.Success"])
        def checkpoint_getMergeConsensus(inputs,outputs):
            output, flagFile = outputs
            runStageCheck('checkpoint', flagFile, outTempPrefix, 'getMergeConsensus')
        stage_count += 1

        if runType == 'pangenome':
            # generate the allele matrix entry for each isolate for each replicon
            def matrixEntryByCoreRep():
                for isolate in full_sequence_list:
                    for repliconName in core_replicons:
                        input = outTempPrefix + 'getRepSNPList/' + refName + '_' + repliconName + '_SNPList.txt'
                        output  = outTempPrefix + isolate + '/deriveRepAlleleMartix/' + refName + '_' + repliconName + '_' + isolate + '_alleles.txt'
                        flagFile = outSuccessPrefix + refName + '_' + repliconName + '_' + isolate + '.deriveRepAlleleMatrix.Success'
                        replicon = repliconName
                        consensus =  outTempPrefix + isolate + '/' + isolate + '_cns.fq'
                        repliconStats = outMerge + refName + '_' + repliconName + '_RepStats.txt'
                        merge_prefix = outMerge
                        yield([input, [output, flagFile], replicon, consensus,repliconStats, merge_prefix])

            @follows(getConsensus, checkpoint_getMergeConsensus)
            @follows(checkpoint_getRepSNPList)
            @files(matrixEntryByCoreRep)
            def deriveRepAlleleMatrix(input, outputs, replicon, consensus, repliconStats, merge_prefix):
                output, flagFile = outputs
                runStageCheck('deriveRepAlleleMatrix', flagFile, input, output, reference, replicon, consensus, repliconStats, merge_prefix)
            stage_count += (len(full_sequence_list)*len(core_replicons))

# checkpoint_deriveRepAlleleMatrix
            @merge(deriveRepAlleleMatrix, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "deriveRepAlleleMatrix.checkpoint.Success"])
            def checkpoint_deriveRepAlleleMatrix(inputs,outputs):
                output, flagFile = outputs
                runStageCheck('checkpoint', flagFile, outTempPrefix, 'deriveRepAlleleMatrix')
            stage_count += 1

            def matrixByCoreRep():
                for repliconName in core_replicons:
                    output  = outTempPrefix + refName + '_' + repliconName + '_alleles.csv'
                    flagFile = outSuccessPrefix + refName + '_' + repliconName + '.collateRepAlleleMatrix.Success'
                    rep_name = refName + '_' + repliconName
                    yield([input==None, [output, flagFile], rep_name])

            @follows(checkpoint_deriveRepAlleleMatrix)
            @files(matrixByCoreRep)
            def collateRepAlleleMatrix(input, outputs, rep_name):
                output, flagFile = outputs
                runStageCheck('collateRepAlleleMatrix', flagFile, outTempPrefix, output, full_sequence_list_string, rep_name)
            stage_count += len(core_replicons)

        else: #runType == 'phylogeny'
            # generate the allele matrix entry for each isolate for each replicon
            def matrixEntryByRep():
                for isolate in full_sequence_list:
                    for repliconName in replicons:
                        input = outTempPrefix + 'getRepSNPList/' + refName + '_' + repliconName[0] + '_SNPList.txt'
                        output  = outTempPrefix + isolate + '/deriveRepAlleleMartix/' + refName + '_' + repliconName[0] + '_' + isolate + '_alleles.txt'
                        flagFile = outSuccessPrefix + refName + '_' + repliconName[0] + '_' + isolate + '.deriveRepAlleleMatrix.Success'
                        replicon = repliconName[0]
                        consensus =  outTempPrefix + isolate + '/' + isolate + '_cns.fq'
                        repliconStats = outMerge + refName + '_' + repliconName[0] + '_RepStats.txt'
                        merge_prefix = outMerge
                        yield([input, [output, flagFile], replicon, consensus,repliconStats, merge_prefix])

            @follows(getConsensus, checkpoint_getMergeConsensus)
            @follows(checkpoint_getRepSNPList)
            @files(matrixEntryByRep)
            def deriveRepAlleleMatrix(input, outputs, replicon, consensus, repliconStats, merge_prefix):
                output, flagFile = outputs
                runStageCheck('deriveRepAlleleMatrix', flagFile, input, output, reference, replicon, consensus, repliconStats, merge_prefix)
            stage_count += (len(full_sequence_list)*len(replicons))

# checkpoint_deriveRepAlleleMatrix
            @merge(deriveRepAlleleMatrix, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "deriveRepAlleleMatrix.checkpoint.Success"])
            def checkpoint_deriveRepAlleleMatrix(inputs,outputs):
                output, flagFile = outputs
                runStageCheck('checkpoint', flagFile, outTempPrefix, 'deriveRepAlleleMatrix')
            stage_count += 1

            def matrixByRep():
                for repliconName in replicons:
                    output  = outTempPrefix + refName + '_' + repliconName[0] + '_alleles.csv'
                    flagFile = outSuccessPrefix + refName + '_' + repliconName[0] + '.collateRepAlleleMatrix.Success'
                    rep_name = refName + '_' + repliconName[0]
                    yield([input==None, [output, flagFile], rep_name])

            @follows(checkpoint_deriveRepAlleleMatrix)
            @files(matrixByRep)
            def collateRepAlleleMatrix(input, outputs, rep_name):
                output, flagFile = outputs
                runStageCheck('collateRepAlleleMatrix', flagFile, outTempPrefix, output, full_sequence_list_string, rep_name)
            stage_count += len(replicons)

        # parse SNP table to create alignment for tree
        @transform(collateRepAlleleMatrix, regex(r"(.*)\/(.+)_alleles.csv"), [outMerge + r"\2_alleles_var_cons"+str(conservation)+".csv", outSuccessPrefix + r"\2_alleles.parseSNPsNoGBK.Success"])
        def parseSNPsNoGBK(inputs, outputs):
            input,_success = inputs
            output, flagFile = outputs
            runStageCheck('parseSNPsNoGBK', flagFile, input, str(conservation), outMerge)
        if runType == "phylogeny":
            stage_count += len(replicons)
        else:
            stage_count += len(core_replicons)

        if conservation != 0.95:
            conservation_temp = 0.95
            @transform(collateRepAlleleMatrix, regex(r"(.*)\/(.+)_alleles.csv"), [outMerge + r"\2_alleles_var_cons" + "0.95.csv", outSuccessPrefix + r"\2_alleles.parseSNPsNoGBK_95.Success"])
            def parseSNPsNoGBK_95(inputs, outputs):
                input,_success = inputs
                output, flagFile = outputs
                runStageCheck('parseSNPsNoGBK', flagFile, input, str(conservation_temp), outMerge)
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)
            
            # create distance matrices based on pair-wise differences in SNPs
            if DifferenceMatrix:
                @follows(parseSNPsNoGBK_95)
                @transform(parseSNPsNoGBK, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outMerge + r"\2_SNP_diff.nxs", outSuccessPrefix + r"\2_alleles.getDifferenceMatrix.Success"])        
                def getDifferenceMatrix(inputs, outputs):
                    output, flagFile = outputs
                    input, _success = inputs
                    runStageCheck('getDifferenceMatrix', flagFile, input)
                if runType == "phylogeny":
                    stage_count += len(replicons)
                else:
                    stage_count += len(core_replicons)

        elif DifferenceMatrix:        
            # create distance matrices based on pair-wise differences in SNPs
            @transform(parseSNPsNoGBK, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outMerge + r"\2_SNP_diff.nxs", outSuccessPrefix + r"\2_alleles.getDifferenceMatrix.Success"])        
            def getDifferenceMatrix(inputs, outputs):
                output, flagFile = outputs
                input, _success = inputs
                runStageCheck('getDifferenceMatrix', flagFile, input)
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)

        # generate tree
        if conservation != 0.95:
            @follows(parseSNPsNoGBK_95)
            @transform(parseSNPsNoGBK, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outMerge + r"\2_alleles_var_cons"+str(conservation)+".tree", outSuccessPrefix + r"\2_alleles.makeTree.Success"])
            def makeTree(inputs, outputs):
                output, flagFile = outputs
                input, _success = inputs
                input = input[:-4] + ".mfasta"
                runStageCheck('makeTree', flagFile, input, output)
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)
        else:
            @transform(parseSNPsNoGBK, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outMerge + r"\2_alleles_var_cons"+str(conservation)+".tree", outSuccessPrefix + r"\2_alleles.makeTree.Success"])
            def makeTree(inputs, outputs):
                output, flagFile = outputs
                input, _success = inputs
                input = input[:-4] + ".mfasta"
                runStageCheck('makeTree', flagFile, input, output)
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)

    else: #refGenbank == False and outMerge == ''

        if runType == 'pangenome':
            # generate the allele matrix entry for each isolate for each replicon
            def matrixEntryByCoreRep():
                for isolate in full_sequence_list:
                    for repliconName in core_replicons:
                        input = outTempPrefix + 'getRepSNPList/' + refName + '_' + repliconName + '_SNPList.txt'
                        output  = outTempPrefix + isolate + '/deriveRepAlleleMartix/' + refName + '_' + repliconName + '_' + isolate + '_alleles.txt'
                        flagFile = outSuccessPrefix + refName + '_' + repliconName + '_' + isolate + '.deriveRepAlleleMatrix.Success'
                        replicon = repliconName
                        consensus =  outTempPrefix + isolate + '/' + isolate + '_cns.fq'
                        repliconStats = outPrefix + refName + '_' + repliconName + '_RepStats.txt'
                        merge_prefix = '-'
                        yield([input, [output, flagFile], replicon, consensus,repliconStats, merge_prefix])

            @follows(getConsensus)
            @follows(checkpoint_getRepSNPList)
            @files(matrixEntryByCoreRep)
            def deriveRepAlleleMatrix(input, outputs, replicon, consensus, repliconStats, merge_prefix):
                output, flagFile = outputs
                runStageCheck('deriveRepAlleleMatrix', flagFile, input, output, reference, replicon, consensus, repliconStats, merge_prefix)
            stage_count += (len(full_sequence_list)*len(core_replicons))

# checkpoint_deriveRepAlleleMatrix
            @merge(deriveRepAlleleMatrix, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "deriveRepAlleleMatrix.checkpoint.Success"])
            def checkpoint_deriveRepAlleleMatrix(inputs,outputs):
                output, flagFile = outputs
                runStageCheck('checkpoint', flagFile, outTempPrefix, 'deriveRepAlleleMatrix')
            stage_count += 1

            def matrixByCoreRep():
                for repliconName in core_replicons:
                    output  = outTempPrefix + refName + '_' + repliconName + '_alleles.csv'
                    flagFile = outSuccessPrefix + refName + '_' + repliconName + '.collateRepAlleleMatrix.Success'
                    rep_name = refName + '_' + repliconName
                    yield([input==None, [output, flagFile], rep_name])

            @follows(checkpoint_deriveRepAlleleMatrix)
            @files(matrixByCoreRep)
            def collateRepAlleleMatrix(input, outputs, rep_name):
                output, flagFile = outputs
                runStageCheck('collateRepAlleleMatrix', flagFile, outTempPrefix, output, full_sequence_list_string, rep_name)
            stage_count += len(core_replicons)

        else: #runType == 'phylogeny'
            # generate the allele matrix entry for each isolate for each replicon
            def matrixEntryByRep():
                for isolate in full_sequence_list:
                    for repliconName in replicons:
                        input = outTempPrefix + 'getRepSNPList/' + refName + '_' + repliconName[0] + '_SNPList.txt'
                        output  = outTempPrefix + isolate + '/deriveRepAlleleMartix/' + refName + '_' + repliconName[0] + '_' + isolate + '_alleles.txt'
                        flagFile = outSuccessPrefix + refName + '_' + repliconName[0] + '_' + isolate + '.deriveRepAlleleMatrix.Success'
                        replicon = repliconName[0]
                        consensus =  outTempPrefix + isolate + '/' + isolate + '_cns.fq'
                        repliconStats = outPrefix + refName + '_' + repliconName[0] + '_RepStats.txt'
                        merge_prefix = '-'
                        yield([input, [output, flagFile], replicon, consensus,repliconStats, merge_prefix])

            @follows(getConsensus)
            @follows(checkpoint_getRepSNPList)
            @files(matrixEntryByRep)
            def deriveRepAlleleMatrix(input, outputs, replicon, consensus, repliconStats, merge_prefix):
                output, flagFile = outputs
                runStageCheck('deriveRepAlleleMatrix', flagFile, input, output, reference, replicon, consensus, repliconStats, merge_prefix)
            stage_count += (len(full_sequence_list)*len(replicons))

# checkpoint_deriveRepAlleleMatrix
            @merge(deriveRepAlleleMatrix, [outTempPrefix+'checkpoint.txt', outSuccessPrefix + "deriveRepAlleleMatrix.checkpoint.Success"])
            def checkpoint_deriveRepAlleleMatrix(inputs,outputs):
                output, flagFile = outputs
                runStageCheck('checkpoint', flagFile, outTempPrefix, 'deriveRepAlleleMatrix')
            stage_count += 1

            def matrixByRep():
                for repliconName in replicons:
                    output  = outTempPrefix + refName + '_' + repliconName[0] + '_alleles.csv'
                    flagFile = outSuccessPrefix + refName + '_' + repliconName[0] + '.collateRepAlleleMatrix.Success'
                    rep_name = refName + '_' + repliconName[0]
                    yield([input==None, [output, flagFile], rep_name])

            @follows(checkpoint_deriveRepAlleleMatrix)
            @files(matrixByRep)
            def collateRepAlleleMatrix(input, outputs, rep_name):
                output, flagFile = outputs
                runStageCheck('collateRepAlleleMatrix', flagFile, outTempPrefix, output, full_sequence_list_string, rep_name)
            stage_count += len(replicons)

        # parse SNP table to create alignment for tree
        @transform(collateRepAlleleMatrix, regex(r"(.*)\/(.+)_alleles.csv"), [outPrefix + r"\2_alleles_var_cons"+str(conservation)+".csv", outSuccessPrefix + r"\2_alleles.parseSNPsNoGBK.Success"])
        def parseSNPsNoGBK(inputs, outputs):
            output, flagFile = outputs
            input,_success = inputs
            runStageCheck('parseSNPsNoGBK', flagFile, input, str(conservation), outPrefix)
        if runType == "phylogeny":
            stage_count += len(replicons)
        else:
            stage_count += len(core_replicons)

        if conservation != 0.95:
            @transform(collateRepAlleleMatrix, regex(r"(.*)\/(.+)_alleles.csv"), [outPrefix + r"\2_alleles_var_cons0.95.csv", outSuccessPrefix + r"\2_alleles.parseSNPsNoGBK_95.Success"])
            def parseSNPsNoGBK_95(inputs, outputs):
                output, flagFile = outputs
                input,_success = inputs
                conservation_temp = 0.95
                runStageCheck('parseSNPsNoGBK', flagFile, input, str(conservation_temp), outPrefix)
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)
            
            # create distance matrices based on pair-wise differences in SNPs
            if DifferenceMatrix:
                @follows(parseSNPsNoGBK_95)
                @transform(parseSNPsNoGBK, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outPrefix + r"\2_SNP_diff.nxs", outSuccessPrefix + r"\2_alleles.getDifferenceMatrix.Success"])        
                def getDifferenceMatrix(inputs, outputs):
                    output, flagFile = outputs
                    input, _success = inputs
                    runStageCheck('getDifferenceMatrix', flagFile, input)
                if runType == "phylogeny":
                    stage_count += len(replicons)
                else:
                    stage_count += len(core_replicons)

        elif DifferenceMatrix:        
            # create distance matrices based on pair-wise differences in SNPs
            @transform(parseSNPsNoGBK, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outPrefix + r"\2_SNP_diff.nxs", outSuccessPrefix + r"\2_alleles.getDifferenceMatrix.Success"])        
            def getDifferenceMatrix(inputs, outputs):
                output, flagFile = outputs
                input, _success = inputs
                runStageCheck('getDifferenceMatrix', flagFile, input)        
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)

        # generate tree
        if conservation != 0.95:
            @follows(parseSNPsNoGBK_95)
            @transform(parseSNPsNoGBK, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outPrefix + r"\2_alleles_var_cons"+str(conservation)+".tree", outSuccessPrefix + r"\2_alleles.makeTree.Success"])
            def makeTree(inputs, outputs):
                output, flagFile = outputs
                input, _success = inputs
                input = input[:-4] + ".mfasta"
                runStageCheck('makeTree', flagFile, input, output)
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)
        else:
            @transform(parseSNPsNoGBK, regex(r"(.*)\/(.+)_alleles_var_cons"+str(conservation)+".csv"), [outPrefix + r"\2_alleles_var_cons"+str(conservation)+".tree", outSuccessPrefix + r"\2_alleles.makeTree.Success"])
            def makeTree(inputs, outputs):
                output, flagFile = outputs
                input, _success = inputs
                input = input[:-4] + ".mfasta"
                runStageCheck('makeTree', flagFile, input, output)
            if runType == "phylogeny":
                stage_count += len(replicons)
            else:
                stage_count += len(core_replicons)

print str(stage_count + 1) + " jobs to be executed in total"
if success_count > 0:
    print str(stage_count+1-success_count) + " jobs left to execute"

# *** Clean up *** 
if outMerge != "":
    if not refGenbank:
        if DifferenceMatrix:
            # delete output directory to finish
            @follows(getDifferenceMatrix, makeTree)
            @files(input==None, outMerge + "finish.deleteDir.Success")
            def deleteDir(input, flagFile):            
                if stage_count > getSuccessCount(outSuccessPrefix):
                    print "\nSome stages have completed without success"
                    print "Pipeline Stopped: check for errors in the log files\n"
                else:
                    make_sequence_list(outMerge, full_sequence_list)
                    make_run_report(outMerge, merge_run, version, run_history, 
                                    pipeline_options.reference, refName, refGenbank, replicons, 
                                    full_sequence_list, readType, runType, core_replicons, 
                                    mapping, bowtie_map_type, replaceReads, minDepth, 
                                    coverFail, depthFail, mappedFail, sdOutgroupMultiplier, 
                                    check_reads_mapped, conservation, modules, bowtie_X_value, 
                                    read_history, sequences)
                    runStageCheck('deleteDir', flagFile, outPrefix)
        else:
            # delete output directory to finish
            @follows(makeTree)
            @files(input==None, outMerge + "finish.deleteDir.Success")
            def deleteDir(input, flagFile):            
                if stage_count > getSuccessCount(outSuccessPrefix):
                    print "\nSome stages have completed without success"
                    print "Pipeline Stopped: check for errors in the log files\n"
                else:
                    make_sequence_list(outMerge, full_sequence_list)
                    make_run_report(outMerge, merge_run, version, run_history, 
                                    pipeline_options.reference, refName, refGenbank, replicons, 
                                    full_sequence_list, readType, runType, core_replicons, 
                                    mapping, bowtie_map_type, replaceReads, minDepth, 
                                    coverFail, depthFail, mappedFail, sdOutgroupMultiplier, 
                                    check_reads_mapped, conservation, modules, bowtie_X_value, 
                                    read_history, sequences)
                    runStageCheck('deleteDir', flagFile, outPrefix)
    else:
        if DifferenceMatrix:
            # delete output directory to finish
            @follows(parseGeneContent, getDifferenceMatrix, makeTree)
            @files(input==None, outMerge + "finish.deleteDir.Success")
            def deleteDir(input, flagFile):
                if stage_count > getSuccessCount(outSuccessPrefix):
                    print "\nSome stages have completed without success"
                    print "Pipeline Stopped: check for errors in the log files\n"
                else:
                    make_sequence_list(outMerge, full_sequence_list)
                    make_run_report(outMerge, merge_run, version, run_history, 
                                    pipeline_options.reference, refName, refGenbank, replicons, 
                                    full_sequence_list, readType, runType, core_replicons, 
                                    mapping, bowtie_map_type, replaceReads, minDepth, 
                                    coverFail, depthFail, mappedFail, sdOutgroupMultiplier, 
                                    check_reads_mapped, conservation, modules, bowtie_X_value, 
                                    read_history, sequences)
                    runStageCheck('deleteDir', flagFile, outPrefix)
        else:
            # delete output directory to finish
            @follows(parseGeneContent, makeTree)
            @files(input==None, outMerge + "finish.deleteDir.Success")
            def deleteDir(input, flagFile):
                if stage_count > getSuccessCount(outSuccessPrefix):
                    print "\nSome stages have completed without success"
                    print "Pipeline Stopped: check for errors in the log files\n"
                else:
                    make_sequence_list(outMerge, full_sequence_list)
                    make_run_report(outMerge, merge_run, version, run_history, 
                                    pipeline_options.reference, refName, refGenbank, replicons, 
                                    full_sequence_list, readType, runType, core_replicons, 
                                    mapping, bowtie_map_type, replaceReads, minDepth, 
                                    coverFail, depthFail, mappedFail, sdOutgroupMultiplier, 
                                    check_reads_mapped, conservation, modules, bowtie_X_value, 
                                    read_history, sequences)
                    runStageCheck('deleteDir', flagFile, outPrefix)    

else:
    if not refGenbank:
        if DifferenceMatrix:
            # delete outTemp directory to finish
            @follows(getDifferenceMatrix, makeTree)
            @files(input==None, outPrefix + "finish.deleteDir.Success")
            def deleteDir(input, flagFile):
                if stage_count > getSuccessCount(outSuccessPrefix):
                    print "\nSome stages have completed without success"
                    print "Pipeline Stopped: check for errors in the log files\n"
                else:
                    make_sequence_list(outPrefix, full_sequence_list)
                    make_run_report(outPrefix, merge_run, version, run_history, 
                                    pipeline_options.reference, refName, refGenbank, replicons, 
                                    full_sequence_list, readType, runType, core_replicons, 
                                    mapping, bowtie_map_type, replaceReads, minDepth, 
                                    coverFail, depthFail, mappedFail, sdOutgroupMultiplier, 
                                    check_reads_mapped, conservation, modules, bowtie_X_value, 
                                    read_history, sequences)
                    runStageCheck('deleteDir', flagFile, outTempPrefix)
        else:
            # delete outTemp directory to finish
            @follows(makeTree)
            @files(input==None, outPrefix + "finish.deleteDir.Success")
            def deleteDir(input, flagFile):
                if stage_count > getSuccessCount(outSuccessPrefix):
                    print "\nSome stages have completed without success"
                    print "Pipeline Stopped: check for errors in the log files\n"
                else:
                    make_sequence_list(outPrefix, full_sequence_list)
                    make_run_report(outPrefix, merge_run, version, run_history, 
                                    pipeline_options.reference, refName, refGenbank, replicons, 
                                    full_sequence_list, readType, runType, core_replicons, 
                                    mapping, bowtie_map_type, replaceReads, minDepth, 
                                    coverFail, depthFail, mappedFail, sdOutgroupMultiplier, 
                                    check_reads_mapped, conservation, modules, bowtie_X_value, 
                                    read_history, sequences)
                    runStageCheck('deleteDir', flagFile, outTempPrefix)
    else:
        if DifferenceMatrix:
            # delete outTemp directory to finish
            @follows(parseGeneContent, getDifferenceMatrix, makeTree)
            @files(input==None, outPrefix + "finish.deleteDir.Success")
            def deleteDir(input, flagFile):
                if stage_count > getSuccessCount(outSuccessPrefix):
                    print "\nSome stages have completed without success"
                    print "Pipeline Stopped: check for errors in the log files\n"
                else:
                    make_sequence_list(outPrefix, full_sequence_list)
                    make_run_report(outPrefix, merge_run, version, run_history, 
                                    pipeline_options.reference, refName, refGenbank, replicons, 
                                    full_sequence_list, readType, runType, core_replicons, 
                                    mapping, bowtie_map_type, replaceReads, minDepth, 
                                    coverFail, depthFail, mappedFail, sdOutgroupMultiplier, 
                                    check_reads_mapped, conservation, modules, bowtie_X_value, 
                                    read_history, sequences)
                    runStageCheck('deleteDir', flagFile, outTempPrefix)
        else:
            # delete outTemp directory to finish
            @follows(parseGeneContent, makeTree)
            @files(input==None, outPrefix + "finish.deleteDir.Success")
            def deleteDir(input, flagFile):
                if stage_count > getSuccessCount(outSuccessPrefix):
                    print "\nSome stages have completed without success"
                    print "Pipeline Stopped: check for errors in the log files\n"
                else:
                    make_sequence_list(outPrefix, full_sequence_list)
                    make_run_report(outPrefix, merge_run, version, run_history, 
                                    pipeline_options.reference, refName, refGenbank, replicons, 
                                    full_sequence_list, readType, runType, core_replicons, 
                                    mapping, bowtie_map_type, replaceReads, minDepth, 
                                    coverFail, depthFail, mappedFail, sdOutgroupMultiplier, 
                                    check_reads_mapped, conservation, modules, bowtie_X_value, 
                                    read_history, sequences)
                    runStageCheck('deleteDir', flagFile, outTempPrefix)

#end of pipeline