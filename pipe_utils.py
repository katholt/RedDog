#!/bin/env python

'''
Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
All rights reserved. (see README.txt for more details)
'''
# Various useful utilities for the pipeline.
# First command (splitPath) is taken from Rubra
import sys
import os.path
import glob
import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def splitPath(path):
    (prefix, base) = os.path.split(path)
    (name, ext) = os.path.splitext(base)
    return (prefix, name, ext)

def isGenbank(refFile):
    ref_file = open(refFile, "rU")
    line = ref_file.readline()
    ref_file.close()
    return line.startswith('LOCUS')   

def isFasta(refFile):
    ref_file = open(refFile, "rU")
    line = ref_file.readline()
    ref_file.close()
    return line.startswith('>')   

def chromInfoFasta(refFile):
    # extract the chromosome name and length from the fasta reference
    chroms = []
    for record in SeqIO.parse(refFile, "fasta"):
        name = record.name
        chroms.append((name, len(record)))
    return chroms

def chromInfoGenbank(refFile):
    # extract the chromosome name and length from the genbank reference
    chroms = []
    for record in SeqIO.parse(refFile, "genbank"):
        name = record.name
        chroms.append((name, len(record)))
    return chroms

# get first value from a simple file
def getValue(coverFile):
    file = open(coverFile, "r")
    values = file.readline()
    file.close()
    valueList = values.split()
    value = valueList[0]
    return value

# get the cover (depth) for a replicon from the coverage file 
# and return as integer * 2
def getCover(coverFile, replicon):
    file = open(coverFile, "r")
    if replicon.find('.') != -1:
        temp_rep = replicon.split('.')
        replicon = temp_rep[0]
    for line in file:
        if line.startswith(replicon):
            values = line.split()
            value = int(float(values[2]) * 2)
    file.close()
    return value

#turn a replicon into a key value(integer) for simple 'hash' table
def get_key(name):
    out_number = 0
    for letter in range(0, len(name)):
        char = name[letter]
        try:
            out = int(char)
        except:
            if (char.upper()=='A') or (char.upper()=='J') or (char.upper()=='S'):
                out = 1
            elif (char.upper()=='B') or (char.upper()=='K') or (char.upper()=='T'):
                out = 2
            elif (char.upper()=='C') or (char.upper()=='L') or (char.upper()=='U'):
                out = 3
            elif (char.upper()=='D') or (char.upper()=='M') or (char.upper()=='V'):
                out = 4
            elif (char.upper()=='E') or (char.upper()=='N') or (char.upper()=='W'):
                out = 5
            elif (char.upper()=='F') or (char.upper()=='O') or (char.upper()=='X'):
                out = 6
            elif (char.upper()=='G') or (char.upper()=='P') or (char.upper()=='Y'):
                out = 7
            elif (char.upper()=='H') or (char.upper()=='Q') or (char.upper()=='Z'):
                out = 8
            elif (char.upper()=='I') or (char.upper()=='R'):
                out = 9
            else:
                out = 9
        out_number += out
    out_number = out_number*(out+1)
    return out_number

# write out the list of sequences to a text file
def make_sequence_list(directory, sequence_list):
    sequence_list_file = open((directory + 'sequence_list.txt') , "w")
    for item in sequence_list:
        sequence_list_file.write(item +'\n')
    sequence_list_file.close()
    return

# count the number of success files in a folder
def getSuccessCount(directory):
    success_files = []
    success_files = glob.glob(directory + '*.Success')
    return len(success_files)

def make_run_report(out_directory,
                    merge_run, 
                    version, 
                    run_history, 
                    reference, 
                    refName, 
                    refGenbank, 
                    replicons, 
                    sequences, 
                    read_type, 
                    run_type, 
                    core_replicons, 
                    mapping, 
                    bowtie_preset, 
                    replace_reads, 
                    min_depth, 
                    coverFail, 
                    depthFail, 
                    mappedFail, 
                    sd_out, 
                    check_reads_mapped, 
                    conservation,
                    modules,
                    bowtie_X,
                    read_history,
                    sequence_list):

    timestamp = str(datetime.datetime.now())
    output = "RedDog " +  version + " " + timestamp +"\n\n"
    output += "Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt\nAll rights reserved. (see README.txt for more details)\n\n"

    if merge_run and run_history != "-":
        output += run_history
        output += "merge\t" + timestamp + "\t" + str(len(sequences)) + "\t" + version + "\n"
    elif merge_run and run_history == "-":
        output += "Run History:\nRun\tdate/time\tsequences\tversion\n"
        output += "merge\t" + timestamp + "\t" + str(len(sequences)) + "\t" + version + "\n"
    else:
        output += "Run History:\nRun\tdate/time\tsequences\tversion\n"
        output += "first\t" + timestamp + "\t" + str(len(sequences)) + "\t" + version + "\n"

    output += "\nModules:\n"
    for module in modules:
        output += "\t" + module + "\n"

    output += "\nReference: " + reference + "\n"
    if refGenbank:
        format = "Genbank"
    else:
        format = "FASTA"
    output += "Reference Format: " + format + "\n"
    output += "No. of Replicons: " + str(len(replicons)) + "\n"
    if run_type == "pangenome":
        output += "Core replicon(s):\n"
        for item in core_replicons:
            output += item + "\n"

    output += "\nNo. of Sequences: " + str(len(sequences)) + "\n"
    if read_type == "PE":
        type_out = "Illumina paired-end"
    elif read_type == "SE":
        type_out = "Illumina single-end"
    else:
        type_out = "Ion Torrent single-end"

    output += "Sequence Type: " + type_out + "\n"
    output += "Sequences List: " + out_directory + "full_sequence_list.txt\n"
    output += "Run type: " + run_type + "\n"

    output += "Mapping: " + mapping + "\n"
    if mapping == 'bowtie':
        output += "bowtie mapping preset: " + bowtie_preset + "\n"
        output += "bowtie X option: " + str(bowtie_X) + "\n"

    if merge_run and replace_reads != '-':
        output += "\nSequences failed by user:\n"
        for name in replace_reads.split(','):
            output += name + "\n"

    output += "\nFilter Options:\n"
    output += "Minimum read depth: " + str(min_depth) + "\n"
    output += "\nPass/fail criteria:\n"
    output += "Coverage of replicon: " + str(coverFail) + "%\n"
    output += "Depth of reads: " + str(depthFail) + "\n"
    if check_reads_mapped != "off":
        output += "Reads mapped: " + str(mappedFail) + "%\n"
        output += "\nReplicon tested for percent of reads mapped:\n"
        output += "Replicon\tPercent of total\n"        
        if len(check_reads_mapped) == 1:
            output += check_reads_mapped[0] + "\t100"
        else:
            found_x = False
            final_ratio = 1.0
            list_of_replicons = []
            ratio_of_replicons = []
            check_reps = check_reads_mapped.split(',')
            for item in check_reps:
                if item != 'x':
                    if found_x != True:
                        list_of_replicons.append(item)
                    else:
                        ratio_of_replicons.append(float(item))
                        final_ratio -= float(item)
                else:
                    found_x = True
            ratio_of_replicons.append(final_ratio)
            for i in range(0, len(list_of_replicons)):
                output += list_of_replicons[i] + "\t" + str(ratio_of_replicons[i]*100) + "\n"
    else:
        output += "Reads mapped: off\n"


    output += "\nAllele conservation ratio: " + str(conservation) + "\n"
    output += "\nOutgroup calling\nStandard Deviations from mean SNP count: " + str(sd_out) + "\n"

    #results by replicon
    if run_type == "pangenome":
        for replicon in core_replicons:
            output += "\nReplicon: " + replicon + "\n"
            stats_file = open((out_directory + refName + '_' + replicon +'_RepStats.txt') , "rU")
            stats_lines = stats_file.readlines()
            failed = 0
            for stats_line in stats_lines:
                if not stats_line.startswith('Isol'):
                    if stats_line.endswith('f\n'):
                        failed += 1
            if failed == 0:
                output += "None of the " +str(len(sequences)) +" isolates failed\n"
            else:
                output += "Of the " +str(len(sequences)) +" isolates, " + str(failed) + " failed\n"

            warnings = []
            warnings = glob.glob(out_directory + replicon + '*_warning.txt')
            if len(warnings) == 1:
                output += "There is one consensus warning file for " + replicon + ":\n"
                output += warnings[0] + "\n"
            elif len(warnings) > 1:
                output += "There are " + str(len(warnings)) + " consensus warning files for " + replicon + ":\n"
                for warning in warnings:
                    output += warning + "\n"

    else:
        for replicon in replicons:
            output += "\nReplicon: " + replicon[0] + "\n"
            stats_file = open((out_directory + refName + '_' + replicon[0] +'_RepStats.txt') , "rU")
            stats_lines = stats_file.readlines()
            failed = 0
            for stats_line in stats_lines:
                if not stats_line.startswith('Isol'):
                    if stats_line.endswith('f\n'):
                        failed += 1
            if failed == 0:
                output += "None of the " + str(len(sequences)) +" isolates failed\n"
            else:
                output += "Of the " + str(len(sequences)) +" isolates, " + str(failed) + " failed\n"

            outgroup_filename = out_directory + refName + '_' + replicon[0] +  '_outgroups.txt'
            if os.path.exists(outgroup_filename):
                outgroup_file = open(outgroup_filename, "rU")
                outgroups = outgroup_file.readlines()
                if len(outgroups) == 1:
                    output += "Outgroup:\n" + outgroups[0]
                else:
                    output += "Outgroups (" +str(len(outgroups)) + "):\n"
                    for outgroup in outgroups:
                        output += outgroup

            warnings = []
            warnings = glob.glob(out_directory + replicon[0] + '_*_warning.txt')
            if len(warnings) == 1:
                output += "There is one consensus warning file for " + replicon[0] + ":\n"
                output += warnings[0] + "\n"
            elif len(warnings) > 1:
                output += "There are " + str(len(warnings)) + " consensus warning files for " + replicon[0] + ":\n"
                for warning in warnings:
                    output += warning + "\n"
    output += "\n"

    if merge_run and read_history != "-":
        output += read_history
        for sequence in sequence_list:
            (prefix, name, ext) = splitPath(sequence)
            output += name + ext + "\t" + prefix + "\t" + timestamp + "\t" + read_type + "\t" + mapping + "\t" + version +"\n"
    elif merge_run and read_history == "-":
        output += "Read History:\nreads\tfull path\tdate/time\tread type\tmapping\tversion\n"
        output += "Note: prior run information unavailable\n"
        for sequence in sequence_list:
            (prefix, name, ext) = splitPath(sequence)
            output += name + ext + "\t" + prefix + "\t" + timestamp + "\t" + read_type + "\t" + mapping + "\t" + version +"\n"
    else:
        output += "Read History:\nreads\tfull path\tdate/time\tread type\tmapping\tversion\n"
        for sequence in sequence_list:
            (prefix, name, ext) = splitPath(sequence)
            output += name + ext + "\t" + prefix + "\t" + timestamp + "\t" + read_type + "\t" + mapping + "\t" + version + "\n"
    output += "\n"

    report_file = open((out_directory + refName + '_run_report.txt') , "w")
    report_file.write(output)
    report_file.close()
    return

def get_run_report_data(run_report):
    report_file = open(run_report, "rU")
    lines = report_file.readlines()
    run_history = ''
    read_history = ''
    ref_name = ''
    ref_format = ''
    replicon_count = ''
    core_replicon_list = []
    run_type = ''
    bowtie_preset = ''
    bowtie_X = ''
    user_failed_list = []
    min_depth = ''
    cover_fail = ''
    depth_fail = ''
    mapped_fail = ''
    replicon_test_list = []
    replicon_percent_list = []
    conservation = ''
    replicon_list = []
    core = False
    read_keep = False
    run_keep = False
    core_keep = False
    user_failed_keep = False
    replicon_list_keep = False

    for line in lines:
        if line == '\n':
            read_keep = False
            run_keep = False
            core_keep = False
            user_failed_keep = False
            replicon_list_keep = False
        if line.startswith('Read History:'):
            read_keep = True
        if read_keep:
            read_history += line
        if line.startswith('Run History:'):
            run_keep = True
        if run_keep:
            run_history += line
        if line.startswith('Reference:'):
            items = line[:-1].split()
            ref_pre, ref_name, ref_ext = splitPath(items[1])
        if line.startswith('Reference Format:'):
            items = line[:-1].split(': ')
            ref_format = items[1]
        if line.startswith('No. of Replicons:'):
            items = line[:-1].split(': ')            
            replicon_count = items[1]
        if line.startswith('Core replicon(s):'):
            core_keep = True
            core = True
        elif core_keep:
            core_replicon_list.append(line[:-1])
        if line.startswith('Run type:'):
            items = line[:-1].split(': ')
            run_type = items[1]
        if line.startswith('bowtie mapping preset:'):
            items = line[:-1].split(': ')
            bowtie_preset = items[1]
        if line.startswith('bowtie X option:'):
            items = line[:-1].split(': ')
            bowtie_X = items[1]
        if line.startswith('Sequences failed by user:'):
            user_failed_keep = True
        elif user_failed_keep:
            user_failed_list.append(line[:-1])
        if line.startswith('Minimum read depth:'):
            items = line[:-1].split(': ')
            min_depth = items[1]
        if line.startswith('Coverage of replicon:'):
            items = line[:-2].split(': ')
            cover_fail = items[1]
        if line.startswith('Depth of reads:'):
            items = line[:-1].split(': ')
            depth_fail = items[1]
        if line.startswith('Reads mapped:'):
            items = line[:-1].split(': ')
            mapped_fail = items[1]
            if mapped_fail.endswith('%'):
                mapped_fail = mapped_fail[:-1]
        if line.startswith('Replicon\tPercent of total'):
            replicon_list_keep = True
        elif replicon_list_keep:
            items = line[:-2].split()
            replicon_test_list.append(items[0])
            replicon_percent_list.append(items[1])
        if line.startswith('Allele conservation ratio:'):
            items = line[:-1].split(': ')
            conservation = items[1]
        if line.startswith('Standard Deviations'):
            items = line[:-1].split(': ')
            sd_out = items[1]
        if not core and line.startswith('Replicon:'):
            items = line[:-1].split()
            replicon_list.append(items[1])

    return (read_history, 
            run_history, 
            ref_name, 
            ref_format, 
            replicon_count, 
            core_replicon_list, 
            run_type, 
            bowtie_preset, 
            bowtie_X, 
            user_failed_list, 
            min_depth, 
            cover_fail, 
            depth_fail, 
            mapped_fail, 
            replicon_test_list, 
            replicon_percent_list, 
            conservation,
            sd_out,
            replicon_list)
 
def getFastaDetails(mfasta):
    mfasta_file = open(mfasta, "rU")
    lines = mfasta_file.readlines()
    isolate_count = 0
    for line in lines:
        if line.startswith('>'):
            isolate_count += 1
    if len(lines) > 1:
        snp_count = len(lines[1])
    else:
        snp_count = 0
    mfasta_file.close()
    return (isolate_count, snp_count)

def checkBases(refFile, fileformat):
    # Read in data
    if fileformat in {'genbank', 'fasta'}:
        seqs = [(record.name, record.seq) for record in SeqIO.parse(refFile, fileformat)]
    else:
        raise ValueError("fileFormat must be either genbank or fasta")
    # Check for nucleotides other than ATGCN
    for desc, seq in seqs:
        nucleotides = set(seq.upper())
        nucleotides_ambiguous = nucleotides - {'A', 'G', 'T', 'C', 'N'}
        if nucleotides_ambiguous:
            sys.stdout.write('error: found ambiguous nucleotides (')
            sys.stdout.write(', '.join(nucleotides_ambiguous))
            sys.stdout.write(') in reference contig ' + desc + '\n')
            sys.exit(1)
