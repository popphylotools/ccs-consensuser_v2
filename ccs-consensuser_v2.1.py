#!/usr/bin/env python3
# coding=utf-8
"""
Run with the -h flag for usage documentation.
"""

import argparse
import copy
import itertools
import logging
import multiprocessing as mp
import subprocess
import os
import glob
import statistics
import sys
from collections import Counter
import csv
import shutil
import re
from zipfile import ZipFile

import Bio
import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# logging configuration
# noinspection PyMissingOrEmptyDocstring
class OneLineExceptionFormatter(logging.Formatter):
    def formatException(self, exc_info):
        result = super().formatException(exc_info)
        return repr(result)

    def format(self, record):
        """Format error message to fit on a single line.
        """
        result = super().format(record)
        if record.exc_text:
            result = result.replace(r"\n", "<newline>")
        return result


handler = logging.StreamHandler()
formatter = OneLineExceptionFormatter(logging.BASIC_FORMAT)
handler.setFormatter(formatter)
root = logging.getLogger()
root.setLevel(os.environ.get("LOGLEVEL", "INFO"))
root.addHandler(handler)

log = logging.getLogger("ccs-consensuser")

# some functions from: https://github.com/chapmanb/bcbb/blob/master/align/adaptor_trim.py
# _remove_adaptor
# trim_adaptor (modified)
# trim_adaptor_w_qual

#Function to check the arguments
def check_args(in_file_list,in_file,oligos,aligner,indel,mode,clean_up,min_prop_reads):
    
    def check_whitespaces(input): 
        if os.path.exists(input):
            with open(input,"r") as file:
                for line in file:
                    if re.search(" ",line.rstrip()):
                        return "whitespace"
        else:
            return "Path"

    if in_file_list is not None:
        result = check_whitespaces(in_file_list)
        if result == "whitespace":
            log.info("Simple spaces are not allowed, please use tab as field separator - {}".format(in_file_list))
            sys.exit(1)
        elif result == "Path":
            log.info("The file does not exist - {}".format(_args.in_file_list))
            sys.exit(1)
    if in_file is not None:
        if re.search("consensus",os.path.basename(in_file)):
            log.info("Please do not use consensus in the fastq file name - {}".format(in_file))
            sys.exit(1)
    if oligos is not None:
        result = check_whitespaces(oligos)
        if result == "whitespace":
            log.info("Simple spaces are not allowed, please use tab as field separator - {}".format(oligos))
            sys.exit(1)
        elif result == "Path":
            log.info("The file does not exist - {}".format(oligos))
            sys.exit(1)
        with open(oligos, 'r') as file:
            reader = csv.reader(file,delimiter="\t")
            for line in reader:
                if line[0] == "forward" or line[0] == "reverse":
                    log.info("This software only supports paired primers format in oligos file - {}".format(oligos))
                    sys.exit(1)

    if aligner not in ["clustalw", "muscle", "mafft"]:
        log.info("Aligner must be clustalw, muscle or mafft")
        sys.exit(1)
    if mode not in ["demultiplex","consensus","consensus_batch","all"]:
        log.info("Mode must be demultiplex, consensus or all")
        sys.exit(1)
    #Check if the external programs are properly installed
    exist = subprocess.call('command -v '+ aligner + '>> /dev/null', shell=True)
    if exist != 0:
        log.info(aligner + "is not installed.")
        sys.exit(1)    
    if mode in ["demultiplex","all"]:
        for program in ["mothur","blastn","makeblastdb"]:
            exist = subprocess.call('command -v '+ program + '>> /dev/null', shell=True)
            if exist != 0:
                log.info(program + "is not installed.")
                sys.exit(1)
    if clean_up not in ["remove","zip","keep"]:
        log.info("clean_up must be remove, zip or keep.")
        sys.exit(1)  
    if indel not in ["NS","N","M"]:
        log.info("indel must be NS, N or M.")
        sys.exit(1)
    if min_prop_reads > 1:
        log.info("min_prop_reads must range between 0 to 1.")
        sys.exit(1)        

    return

#Convert a relative path to full path
def convert_to_full_path(path):
    absolute_path = os.path.abspath(path)
    return absolute_path

#It belongs to demultiplex_wraper function
#Extract the primers from the oligos file (mothur) and returns a list of primerF and primerR, a list of sample IDs, primer names (if more than one marker) as well as a fasta file.
def parse_oligos(oligos,output_dir):
    
    primer_fasta_path = os.path.join(output_dir,"primers.fasta")
    primer_info = []
    primer_names = []
    barcode_info = []
    with open(oligos, 'r') as file:
        reader = csv.reader(file,delimiter="\t")
        for line in reader:
            if line[0] == "primer" and len(line) == 3:
                primer_info.append([line[1],line[2]])
                primer_names.append(["markerID",line[1],line[2]])
            elif line[0] == "primer" and len(line) == 4:
                primer_info.append([line[1],line[2]])
                primer_names.append([line[3],line[1],line[2]])
            elif line[0] == "barcode" and len(line) == 4:
                barcode_info.append(line[3])

    if len(primer_info) == 0 or len(barcode_info) == 0:
        log.info("Check your oligos file - {}".format(oligos))
        return
    if len(primer_names) != len(primer_info):
        log.info("Check your oligos file. You may include primers ID - {}".format(oligos))
        return
    #Create ids
    count = 1
    primer_info_new = []
    primerF = []
    primerR = []
    for primer in primer_info:
        #For primer forward
        #Create a list with all possible primers based on the ambiguous nucleotides
        list_primers = primers_ambiguity(primer[0])
        counter_ambiguity_F = 1
        for primer_no_ambiguity in list_primers:
            barcode_name_F = "primerF" + "_" + str(count) + "_" + str(counter_ambiguity_F)
            primer_info_new.append([barcode_name_F,primer_no_ambiguity])
            primerF.append(barcode_name_F)
            counter_ambiguity_F = counter_ambiguity_F + 1

        #For primer reverse
        #Create a list with all possible primers based on the ambiguous nucleotides
        list_primers = primers_ambiguity(primer[1])
        counter_ambiguity_R = 1
        for primer_no_ambiguity in list_primers:
            barcode_name_R = "primerR" + "_" + str(count) + "_" + str(counter_ambiguity_R)
            primer_info_new.append([barcode_name_R,primer_no_ambiguity])
            primerR.append(barcode_name_R)
            counter_ambiguity_R = counter_ambiguity_R + 1
        count += 1

    with open(primer_fasta_path, 'w') as primer_fasta:
        for primer in primer_info_new:
            primer_fasta.write(">" + primer[0] + "\n" + primer[1] + "\n")
    return primerF,primerR,barcode_info,primer_names,primer_fasta_path

#It belongs to demultiplex_wraper function
#This function is part of parse_oligos and primers_ambiguity functions
#Creates all possible primer combinations determined by the ambiguous bases on the primer
def split_primers(primer):
    dict_ambiguity = {"A": ["A"],"C":["C"],"G":["G"],"T":["T"],"R":["A","G"],"Y":["C","T"],"S":["G","C"],"W":["A","T"],"K":["G","T"],"M":["A","C"],"B":["C","G","T"],"D":["A","G","T"],"H":["A","C","T"],"V":["A","C","G"],"N":["A","T","C","G"],"-":["-"]}

    amb_index = []
    amb_nuc = []
    new_primers = []
    primer_list = list(primer)
    for idx in range(len(primer)):
        if primer[idx] in dict_ambiguity.keys():
            amb_index.append(idx)
            amb_nuc.append(dict_ambiguity[primer[idx]])
    combinations = list(itertools.product(*amb_nuc))
    for comb in combinations:
        current_primer_list = copy.deepcopy(primer_list)
        for comb_amb_index in range(len(comb)):
            current_comb_amb_index = amb_index[comb_amb_index]
            current_primer_list[current_comb_amb_index] = comb[comb_amb_index]
        new_primers.append(''.join(current_primer_list))
    return new_primers

#It belongs to demultiplex_wraper function
#This funxtion is part of parse_oligos function
#Test if there is an ambiguous bases to invoke split_primers function.
def primers_ambiguity(primer):
    Ambiguous_bases=["R","Y","S","W","K","M","B","D","H","V","N"]
    if len(primer) > 20:
        primer = primer[:20]
    for amb_base in Ambiguous_bases:
        if amb_base in primer:
            list_new_primers = split_primers(primer)
            return list_new_primers
    list_primers = [primer]
    return list_primers

#It belongs to demultiplex_wraper function
#Convert a fastq to fasta file. It returns the fasta path and number of reads.
def convert_fastq2fasta(input_fastq,output_dir):
    output_basename_fasta = os.path.splitext(os.path.basename(input_fastq))[0] + ".fasta"
    output_fasta = os.path.join(output_dir,output_basename_fasta)
    with open(input_fastq, "r") as input_handle:
        with open(output_fasta, "w") as output_handle:
            count = 0
            for record in SeqIO.parse(input_handle, "fastq"):
                SeqIO.write(record, output_handle, "fasta")
                count = count+1
            number_seq = count

    return output_fasta,number_seq

#It belongs to demultiplex_wraper function
#Align the primers to the reads using blast-short tool. It returns the path of the blast output.
def short_blast(input_fasta,output_dir,number_seq,primer_fasta,number_processors):
    output_blast = os.path.join(output_dir,"blastn_short_primer_vs_input_fasta")
    blast_db = os.path.join(output_dir,"reads_blast_db")
    try:
        cmd = ["makeblastdb", "-in", input_fasta, "-input_type","fasta","-dbtype","nucl","-out",blast_db]
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        p.wait()

    except OSError as e:
        if e.errno == errno.ENOENT:
        # handle file not found error.
            log.info("There was an error while running makeblastdb command - {}".format(e.errno))
        else:
            # Something else went wrong while trying to run "makeblastdb"
            raise
    try:
        cmd = ["blastn","-task","blastn-short","-db",blast_db,"-query",primer_fasta,"-out",output_blast,"-evalue","10","-outfmt","6","-max_target_seqs",str(number_seq*5),"-num_threads",str(number_processors)]
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        p.wait()

    except OSError as e:
        if e.errno == errno.ENOENT:
        # handle file not found error.
            log.info("There was an error while running blastn command  - {}".format(e.errno))
        else:
        # Something else went wrong while trying to run "blastn"
            raise

    #Clean the temporary files
    try:
        os.remove(input_fasta)
        output_consensus_files = glob.glob(blast_db + "*")
        for db_intermediate in output_consensus_files:
            os.remove(db_intermediate)
    
    except OSError as e:  ## if failed, report it back to the user ##
        log.info("There was an error removing the blastn intermidiate files, {} - {}".format(e.filename, e.strerror))

    return output_blast

#It belongs to demultiplex_wraper function
#Parse the blast results to generate a dictionary (key = read_id, value = orientation).
def convert_blastoutfmt6_2_dict(output_blast,list_primerF,list_primerR):
    blast_result_dict = {}
    with open(output_blast, 'r') as file:
        reader = csv.reader(file,delimiter="\t")
        for row in reader:
            if row[1] not in blast_result_dict.keys():
                list_info = [row[0], int(row[6]),int(row[7]),int(row[8]),int(row[9])]
                orientation = test_orientation(list_info,list_primerF,list_primerR)
                if orientation is not None:
                    blast_result_dict[row[1]] = orientation
                else:
                    log.info("Something went wrong in convert_blastoutfmt6_2_dict function trying to get the orientation of the reads")
    #Clean the temporary files
    try:
        os.remove(output_blast)

    except OSError as e:  ## if failed, report it back to the user ##
        log.info("There was an error removing the blastn result, {} - {}".format(output_blast, e.strerror))

    return blast_result_dict

#It belongs to demultiplex_wraper function
#Determine the orientation of the read based on a list that contains the query_ID (column 1) and coordinates (columns 7, 8, 9 and 10) etracted from the blast output.
def test_orientation(list,list_primerF, list_primerR):
    if list[0] in list_primerF:
        if list[2] > list[1] and list[4] > list[3]:
            return "oriented_F"
        elif list[2] < list[1] and list[4] < list[3]:
            return "oriented_F"
        elif list[2] > list[1] and list[4] < list[3]:
            return "oriented_R"
        elif list[2] < list[1] and list[4] > list[3]:
            return "oriented_R"
    elif list[0] in list_primerR:
        if list[2] > list[1] and list[4] > list[3]:
            return "oriented_R"
        elif list[2] < list[1] and list[4] < list[3]:
            return "oriented_R"
        elif list[2] > list[1] and list[4] < list[3]:
            return "oriented_F"
        elif list[2] < list[1] and list[4] > list[3]:
            return "oriented_F"

#It belongs to demultiplex_wraper function
#Create a fastq file oriented based on the primerF.
def reorientation(input_fastq,output_dir,orientation_dict):
    output_oriented_basename_fastq = os.path.splitext(os.path.basename(input_fastq))[0] + "_oriented.fastq"
    output_oriented_fastq = os.path.join(output_dir,output_oriented_basename_fastq)
    
    with open(output_oriented_fastq, "w") as oriented_fastq:
        for record in SeqIO.parse(input_fastq,'fastq'):
            if record.id in orientation_dict.keys():
                orientation = orientation_dict[record.id]
                if orientation == "oriented_F":
                    SeqIO.write(record, oriented_fastq, 'fastq')
                elif orientation == "oriented_R":
                    record_rc = record.reverse_complement()
                    record_rc.id = record.id
                    record_rc.description = ''
                    SeqIO.write(record_rc, oriented_fastq, 'fastq')
                else:
                    log.info("There was an error while running the reorientation function")

    return output_oriented_fastq

#It belongs to demultiplex_wraper function
#This function runs the fastq.info of the mothur program to demultiplex the reads.
def run_mothur(oriented_fastq,oligos_path,primer_mismatch,barcode_mismatch):
    output_dir = os.path.dirname(oriented_fastq) 
    command = "#fastq.info(fastq=" + oriented_fastq + ",oligos=" + oligos_path + ",pacbio=T,qfile=F,fasta=F,bdiffs=" + str(barcode_mismatch) + ",pdiffs=" + str(primer_mismatch) + ")"
    
    try:
        cmd = ["mothur",command]
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,cwd=output_dir)
        p.wait()

    except OSError as e:
        if e.errno == errno.ENOENT:
        # handle file not found error.
            log.info("There was an error while running the mothur command - {}".format("mothur "+command))
        else:
        # Something else went wrong while trying to run "mothur"
            raise
            return
    return

#This function accomodates the orientation of the reads and runs the fastq.info of the mothur program to demultiplex the reads.
def demultiplex_wraper(input_fastq,oligos_path,output_dir,number_processors,primer_mismatch,barcode_mismatch):
    list_primerF,list_primerR,barcode_info,primer_names,primer_fasta_path = parse_oligos(oligos_path,output_dir)
    input_fasta,number_seq = convert_fastq2fasta(input_fastq,output_dir)
    output_blast = short_blast(input_fasta,output_dir,number_seq,primer_fasta_path,number_processors)
    orientation_dict = convert_blastoutfmt6_2_dict(output_blast,list_primerF,list_primerR)
    output_oriented_fastq = reorientation(input_fastq,output_dir,orientation_dict)
    run_mothur(output_oriented_fastq,oligos_path,primer_mismatch,barcode_mismatch)
    return output_oriented_fastq,barcode_info,primer_names

#Belongs to trim_and_mask_seq_records.
def _remove_adaptor(seq, region, right_side=True):
    """Remove an adaptor region and all sequence to the right or left.
    """
    if right_side:
        try:
            pos = seq.find(region)
        # handle Biopython SeqRecords
        except AttributeError:
            pos = seq.seq.find(region)
        return seq[:pos]
    else:
        try:
            pos = seq.rfind(region)
        # handle Biopython SeqRecords
        except AttributeError:
            pos = seq.seq.rfind(region)
        return seq[pos + len(region):]

#Belongs to trim_and_mask_seq_records.
def trim_adaptor(seq, adaptor, primer_mismatch, right_side=True):
    """Trim the given adaptor sequence from a starting sequence.
    * seq can be either of:
       - string
       - Seq
    * adaptor is a string sequence
    * primer_mismatch specifies how many errors are allowed in the match between
    adaptor and the base sequence. Matches with more than this number of errors
    are not allowed.
    """
    gap_char = '-'
    seq_original = seq
    seq = str(seq).upper()
    adaptor = str(adaptor).upper()
    exact_pos = seq.find(adaptor)
    if exact_pos >= 0:
        seq_region = str(seq[exact_pos:exact_pos + len(adaptor)])
        adapt_region = adaptor
    else:
        aligner_local = PairwiseAligner(mode = 'local', match_score = 5.0, mismatch_score = -4.0, open_gap_score = -9.0, extend_gap_score = -0.5)
        aligns = aligner_local.align(seq, adaptor)
        if len(aligns) == 0:
            adapt_region, seq_region = ("", "")
        else:
            seq_region, adapt_region = aligns[0]
    matches = match_func_amb(seq_region,adapt_region)
    # too many errors -- no trimming
    if (len(adaptor) - matches) > primer_mismatch:
        return seq_original
    # remove the adaptor sequence and return the result
    else:
        return _remove_adaptor(seq_original, seq_region.replace(gap_char, ""),right_side)

#Belongs to trim_and_mask_seq_records.
#Determine the number of matches taking into account ambiguous bases.
def match_func_amb(seq_region,adapt_region):
    dict_ambiguity = {"A": ["A"],"C":["C"],"G":["G"],"T":["T"],"R":["A","G"],"Y":["C","T"],"S":["G","C"],"W":["A","T"],"K":["G","T"],"M":["A","C"],"B":["C","G","T"],"D":["A","G","T"],"H":["A","C","T"],"V":["A","C","G"],"N":["A","T","C","G"],"-":["-"]}
    count = 0
    for idx,s in enumerate(seq_region):
        base_adapt = adapt_region[idx]
        possible_bases_adapt = dict_ambiguity[base_adapt]
        if s in possible_bases_adapt:
            count = count + 1
    return count


#Belongs to trim_and_mask_seq_records.
def trim_adaptor_w_qual(seq, qual, adaptor, primer_mismatch, right_side=True):
    """Trim an adaptor with an associated quality string.
    Works like trimmed adaptor, but also trims an associated quality score.
    """
    assert len(seq) == len(qual)
    tseq = trim_adaptor(seq, adaptor, primer_mismatch, right_side=right_side)
    if right_side:
        pos = seq.find(tseq)
    else:
        pos = seq.rfind(tseq)
    tqual = qual[pos:pos + len(tseq)]
    assert len(tseq) == len(tqual)
    return tseq, tqual

def gap_consensus(input_consensus, threshold=.7, mask_char="N", consensus_ambiguous_character="n",
                  require_multiple=False, consensus_ignore_mask_char=False):
    """Output a fast consensus sequence of the alignment, allowing gaps.
    Same as dumb_consensus(), but allows gap on the output.
    Things to do:
     - Let the user define that with only one gap, the result
       character in consensus is gap.
     - Let the user select gap character, now
       it takes the same as input.
    """
    # Iddo Friedberg, 1-JUL-2004: changed ambiguous default to "X"
    
    # parse aligned fasta
    with open(input_consensus, "r") as f:
        alignment = AlignIO.read(f, "fasta")

    # take consensus of aligned fasta
    summary_align = AlignInfo.SummaryInfo(alignment)

    consensus = ''

    # find the length of the consensus we are creating
    con_len = summary_align.alignment.get_alignment_length()

    # go through each seq item
    for n in range(con_len):
        # keep track of the counts of the different atoms we get
        atom_dict = {}
        num_atoms = 0

        for record in summary_align.alignment:
            # make sure we haven't run past the end of any sequences
            # if they are of different lengths
            if n < len(record.seq):
                if consensus_ignore_mask_char and (record.seq[n] == mask_char):
                    continue
                if record.seq[n] not in atom_dict:
                    atom_dict[record.seq[n]] = 1
                else:
                    atom_dict[record.seq[n]] += 1

                num_atoms += 1

        max_atoms = []
        max_size = 0

        for atom in atom_dict:
            if atom_dict[atom] > max_size:
                max_atoms = [atom]
                max_size = atom_dict[atom]
            elif atom_dict[atom] == max_size:
                max_atoms.append(atom)

        if require_multiple and num_atoms == 1:
            consensus += consensus_ambiguous_character
        elif (len(max_atoms) == 1) and ((float(max_size) / float(num_atoms)) >= threshold):
            consensus += max_atoms[0]
        else:
            consensus += consensus_ambiguous_character
    return Seq(consensus)

def create_unique_dir(path, limit=99):
    """Return a path to an empty directory. Either the dir at path, or a dir of the form 'path + _01'
    :param path: The initial path to use
    :param limit: The maximum number of directory variations this function will attempt to create.
    :return: A path to an empty directory.
    """
    width = len(str(limit))
    original = path.rstrip(os.sep)
    if len(os.listdir(original)) == 0:
        log.info("Using output directory - {}".format(path))
        return original  # folder empty, let's use it
    count = 1
    while count < limit:
        try:
            os.mkdir(path)
            log.info("Creating output directory - {}".format(path))
            return path
        except OSError as path_error:
            if path_error.errno == 17:  # file exists
                path = "{0}_{1:0>{2}}".format(original, count, width)
                count += 1
            else:
                raise
    else:
        msg = "could not uniquely create directory {0}: limit `{1}` reached"
        raise Exception(msg.format(original, limit))


#Belongs to trim_and_mask_seq_records.
def trim_both_ends(seq_rec, primer_a, primer_b, primer_mismatch, reverse_complement=False):
    if reverse_complement:
        rc = seq_rec.reverse_complement()
        un_trimed_seq = rc.seq
        un_trimed_qual = rc.letter_annotations["phred_quality"]
    else:
        un_trimed_seq = seq_rec.seq
        un_trimed_qual = seq_rec.letter_annotations["phred_quality"]

    # primer A,B found
    found_a = False
    found_b = False

    half_trimed_seq, half_trimed_qual = trim_adaptor_w_qual(un_trimed_seq,
                                                            un_trimed_qual,
                                                            adaptor=primer_a,
                                                            primer_mismatch=primer_mismatch,
                                                            right_side=False)
    if len(half_trimed_seq) < len(un_trimed_seq):
        found_a = True

    full_trimed_seq, full_trimed_qual = trim_adaptor_w_qual(half_trimed_seq,
                                                            half_trimed_qual,
                                                            adaptor=primer_b,
                                                            primer_mismatch=primer_mismatch,
                                                            right_side=True)
    if len(full_trimed_seq) < len(half_trimed_seq):
        found_b = True

    if found_a and found_b:
        trimed_seq_rec = copy.deepcopy(seq_rec)
        del trimed_seq_rec.letter_annotations["phred_quality"]
        trimed_seq_rec.seq = full_trimed_seq
        trimed_seq_rec.letter_annotations["phred_quality"] = full_trimed_qual
        return trimed_seq_rec
    else:
        return None

#Belongs to trim_and_mask_seq_records.
def mask_seq_record(seq_rec, min_score, mask_char, inplace=False):
    if not inplace:
        masked_seq_req = copy.deepcopy(seq_rec)
    else:
        masked_seq_req = seq_rec

    base_list = list(masked_seq_req.seq)
    for loc in range(len(base_list)):
        if masked_seq_req.letter_annotations["phred_quality"][loc] < min_score:
            base_list[loc] = mask_char

    masked_seq_req.seq = Seq("".join(base_list))

    return masked_seq_req

#Function to trim and mask bases based on phred quality
def trim_and_mask_seq_records(records, primer_a, primer_b, primer_mismatch, min_base_score, basename, mask_char,
                              min_seq_score=None):
    for seq_rec in records:
        trimed_seq_rec = trim_both_ends(seq_rec, primer_a, primer_b, primer_mismatch, reverse_complement=False)

        # primers found in forward direction
        if trimed_seq_rec is not None and len(trimed_seq_rec.letter_annotations["phred_quality"]) >= 2:
            avg_score = np.mean(trimed_seq_rec.letter_annotations["phred_quality"])
            if min_seq_score and (avg_score < min_seq_score):
                log.info("seq excluded - avg_score:{:4.2f} < min_seq_score:{} - {} {}".format(avg_score, min_seq_score,
                                                                                              basename, seq_rec.id))
                continue
            yield mask_seq_record(trimed_seq_rec, min_base_score, mask_char=mask_char, inplace=True)

        # primers not found in forward direction
        else:
            trimed_seq_rec = trim_both_ends(seq_rec, primer_a, primer_b, primer_mismatch, reverse_complement=True)
            if trimed_seq_rec is not None and len(trimed_seq_rec.letter_annotations["phred_quality"]) >= 2:  # primers found in reverse direction
                avg_score = np.mean(trimed_seq_rec.letter_annotations["phred_quality"])
                if min_seq_score and (avg_score < min_seq_score):
                    log.info(
                        "seq excluded - avg_score:{:4.2f} < min_seq_score:{} - {} {}".format(avg_score, min_seq_score,
                                                                                             basename, seq_rec.id))
                    continue
                yield mask_seq_record(trimed_seq_rec, min_base_score, mask_char=mask_char, inplace=True)

            # primers not found in either direction
            else:
                log.info("seq excluded - primers not found - {} {}".format(basename, seq_rec.id))
                continue

#Function to produce multiple sequence alignment
def alignment_step(fasta_path,output_path,aligner,number_seqs,basename):
    
    if aligner == "muscle":
        output_dir = os.path.dirname(output_path)
        infile = fasta_path
        outfile = output_path
        try:
            cmd = ["muscle", "-in", infile,"-out",outfile,"-quiet"]
            p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL,stderr=subprocess.PIPE, cwd=output_dir)
            p.wait()
            (stdout, stderr) = p.communicate()
        except calledProcessError as err:
            log.info("alignment failed - {} - {}".format(err.stderr, basename))
            return 

    elif aligner == "clustalw":
        output_dir = os.path.dirname(output_path)
        infile = fasta_path
        outfile_nexus = os.path.splitext(output_path)[0] + ".nex"
        try:
            cmd = ["clustalw2", "-INFILE=" + infile,"-TYPE=DNA","-OUTPUT=NEXUS","-OUTFILE=" + outfile_nexus]
            p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL,stderr=subprocess.PIPE, cwd=output_dir)
            p.wait()
            (stdout, stderr) = p.communicate()
        except calledProcessError as err:
            log.info("alignment failed - {} - {}".format(err.stderr, basename))
            return        
        
        #Convert nexus to fasta
        outfile_fasta = os.path.splitext(outfile_nexus)[0] + ".fasta"
        with open(outfile_nexus, "r") as input_handle:
            with open(outfile_fasta, "w") as output_handle:
                sequences = SeqIO.parse(input_handle, "nexus")
                count = SeqIO.write(sequences, output_handle, "fasta")

        file = os.path.splitext(infile)[0] + ".dnd"
        try:
            os.remove(file)
            os.remove(outfile_nexus)
        except:
            log.info("Error while deleting clustal tree file in alignment_step function: {}".format(file))

    elif aligner == "mafft":
        output_dir = os.path.dirname(output_path)
        infile = fasta_path
        outfile = output_path
        if number_seqs < 200:
            with open(outfile,"wb") as out_handle:
                try:
                    cmd = ["mafft", "--localpair","--maxiterate","1000","--thread","1","--quiet","--preservecase", infile]
                    p = subprocess.Popen(cmd, stdout=out_handle,stderr=subprocess.PIPE, cwd=output_dir)
                    p.wait()
                    (stdout, stderr) = p.communicate()
                except calledProcessError as err:
                    log.info("alignment failed - {} - {}".format(err.stderr, basename))
                    return
        
        else:
            with open(outfile,"wb") as out_handle:
                try:
                    cmd = ["mafft", "--auto","--thread","1","--quiet","--preservecase", infile]
                    p = subprocess.Popen(cmd, stdout=out_handle,stderr=subprocess.PIPE, cwd=output_dir)
                    p.wait()
                    (stdout, stderr) = p.communicate()
                except calledProcessError as err:
                    log.info("alignment failed - {} - {}".format(err.stderr, basename))
                    return
    
    else:
        log.info("alignment failed - invalid aligner: {} - {}".format(aligner, basename))
        return

#Main function to mask bases based on coverage (minimum count). Returns Biopython alingment object.
def mask_min_count_wraper(aln_fasta,min_SNP_cov,mask_char,indel):
    
    aln_dna = AlignIO.read(open(aln_fasta), 'fasta')

    dict_pos={}
    for python_pos in range(aln_dna.get_alignment_length()):
        pos_aln = aln_dna[:,python_pos].upper()
        base_counter = Counter(pos_aln)
        low_coverage = analize_counter(base_counter,min_SNP_cov,indel)
        if len(low_coverage) > 0:
            dict_pos[python_pos] = low_coverage

    records = mask_low_coverage(aln_dna,mask_char,dict_pos)

    return records

#Return the characters with low coverage (based on min_SNP_cov). Belongs to mask_min_count_wraper(aln_fasta,min_SNP_cov).
def analize_counter(base_counter,min_SNP_cov,indel):
    low_coverage = []
    bases = ["A","T","C","G"]
    if indel == "T":
        bases.append("-")
    for base in base_counter.keys():
        if base_counter[base] < min_SNP_cov and base in bases:
            low_coverage.append(base)

    return low_coverage

#Mask the bases based on a dictionary (key = python position and values list of bases to be masked). Belongs to mask_min_count_wraper(aln_fasta,min_SNP_cov).
def mask_low_coverage(aln_dna,mask_char,dict_pos):
    for idx_record in range(len(aln_dna)):
        seq = str(aln_dna[idx_record].seq)
        
        new_seq = seq[:]
        for base_pos in range(len(seq) - 1, -1, -1):
            if base_pos in dict_pos.keys() and seq[base_pos] in dict_pos[base_pos]:
                new_seq = new_seq[:base_pos] + mask_char + new_seq[base_pos+1:]
        aln_dna[idx_record].seq = Seq(new_seq)

    return aln_dna

#Wraper to get the "haplotypes" (unique sequences). It returns a list of path to fasta unaligned sequences and haplotype frequency (number of seqs per file).
def wraper_get_unique_seqs(path_fasta,basename,aligner,min_seq_cov,indel):
    records = list(SeqIO.parse(path_fasta, "fasta"))
    
    list_record_ids = [record.id for record in records]
    dict_unique ={}
    for record1,record2  in itertools.combinations(records, 2):
        if record1.id not in dict_unique.keys():
            dict_unique[record1.id]=[record1.id]
        if record2.id != None:
            seq1 = str(record1.seq)
            seq2 = str(record2.seq)
            matches = compare_seq(seq1,seq2,indel)
            if matches == len(seq1):
                dict_unique[record1.id].append(record2.id)
                index = list_record_ids.index(record2.id)
                records[index].id = None
    
    #Clean the dictionary
    if None in dict_unique.keys():
        del dict_unique[None]

    new_fasta_path_freq_list = dict_unique_to_files(path_fasta, basename, dict_unique,min_seq_cov)

    return new_fasta_path_freq_list

#This function belongs to wraper_get_unique_seqs
#Compare two sequences and return the number of matches between sequences. Missing data (N) is considered a wildcard, so it will match any base or even a gap.
def compare_seq(seq1,seq2,indel):
    counter = 0
    #if NS, gaps are treated as a new state. This may be useful to resolve indels (variants), but may require more sequencing coverage to resolve the haplotypes (old F). 
    if indel == "NS":
        for base_idx in range(len(seq1)):
            if seq1[base_idx] == seq2[base_idx]:
                counter = counter + 1
                continue
            elif (seq1[base_idx] == "N" and seq2[base_idx] == "-") or (seq1[base_idx] == "-" and seq2[base_idx] == "N"):
                counter = counter
                continue
            elif seq1[base_idx] == "N" or seq2[base_idx] == "N":
                counter = counter + 1
                continue
    
    #If N, gaps are treated as a fifth nucleotide (old T). 
    elif indel == "N":
        for base_idx in range(len(seq1)):
            if seq1[base_idx] == seq2[base_idx]:
                counter = counter + 1
                continue
            elif seq1[base_idx] == "N" or seq2[base_idx] == "N":
                counter = counter + 1
                continue

    #If M, gaps are treated as missing data (ignored).
    elif indel == "M":
        for base_idx in range(len(seq1)):
            if seq1[base_idx] == seq2[base_idx]:
                counter = counter + 1
                continue
            elif seq1[base_idx] == "N" or seq2[base_idx] == "N":
                counter = counter + 1
                continue
            elif seq1[base_idx] == "-" or seq2[base_idx] == "-":
                counter = counter + 1
                continue
    return counter

#This function belongs to wraper_get_unique_seqs
#Convert a dictionary with haplotypes (id1:[id1,id2,...]) to separate files.
def dict_unique_to_files(path_fasta_aln, basename, dict_unique,min_seq_cov):
    output_dir = os.path.dirname(path_fasta_aln)
    records_dict = SeqIO.to_dict(SeqIO.parse(path_fasta_aln, "fasta"))
    counter = 1
    new_fasta_path_freq_list = []
    for haplotype in dict_unique.keys():
        haplotype_records = []
        for similar_seqs_id in dict_unique[haplotype]:
            #Populate the temporal list of records
            record = records_dict[similar_seqs_id]
            record.seq = Seq(remove_gaps(str(record.seq)))
            haplotype_records.append(record)
        
        number_seq = len(haplotype_records)
        #Create a path
        new_fasta_path = os.path.join(output_dir, basename) + "_" + str(counter) + "." + str(number_seq) + ".fasta"

        with open(new_fasta_path, "wt") as f:
            for r in haplotype_records:
                f.write(r.format("fasta"))
        counter = counter + 1
        
        new_fasta_path_freq_list.append([new_fasta_path,number_seq])

    return new_fasta_path_freq_list

#Create a directory
def create_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    return dir_path

#It belongs to all mode.

def prepare_directories_files_all_mode(in_file,output_dir,barcode_info,primer_names,oligos):
    
    #For more multiple primer pairs (markers)
    input_paths = []
    primerF = []
    primerR = []
    new_output_dirs = []
    dict_markers_outputDIR = {}
    if primer_names[0][0] != "markerID":
        #Create a directory for each marker
        for marker_info in primer_names:
            marker_id = marker_info[0]
            marker_dir = os.path.join(output_dir,marker_id)
            create_dir(marker_dir)
        
        oriented_basename = os.path.splitext(os.path.basename(in_file))[0] + "_oriented"
        for sample_marker in itertools.product(barcode_info, primer_names):
            sample_id = sample_marker[0]
            marker_id = sample_marker[1][0]
            marker_info_primerF = sample_marker[1][1]
            marker_info_primerR = sample_marker[1][2]
            output_oriented_demultiplex_fastq = os.path.join(output_dir,oriented_basename) + "." + sample_id + "." + marker_id + ".fastq"
            if os.path.exists(output_oriented_demultiplex_fastq):
                output_sample_dir = os.path.join(output_dir,marker_id,sample_id)
                create_dir(output_sample_dir)
                new_fastq_sample_path = os.path.join(output_sample_dir,oriented_basename) + "." + sample_id + "." + marker_id + ".fastq"
                shutil.move(output_oriented_demultiplex_fastq, new_fastq_sample_path)
                input_paths.append(new_fastq_sample_path)
                primerF.append(marker_info_primerF)
                primerR.append(marker_info_primerR)
                new_output_dirs.append(output_sample_dir)
                #Populate the dict_markers_outputDIR
                if marker_id not in dict_markers_outputDIR.keys():
                    dict_markers_outputDIR[marker_id] = []
                    dict_markers_outputDIR[marker_id].append(output_sample_dir)
                elif marker_id in dict_markers_outputDIR.keys():
                    dict_markers_outputDIR[marker_id].append(output_sample_dir)

        #Remove empty marker folders
        for marker_info in primer_names:
            marker_directory = os.path.join(output_dir,marker_info[0])
            if len(os.listdir(marker_directory)) == 0:
                os.rmdir(marker_directory)    
        return input_paths,primerF,primerR,new_output_dirs,dict_markers_outputDIR
    
    #For one primer pair (markers)
    elif len(primer_names) == 1 and primer_names[0][0] == "markerID":
        oriented_basename = os.path.splitext(os.path.basename(in_file))[0] + "_oriented"
        marker_info_primerF = primer_names[0][1]
        marker_info_primerR = primer_names[0][2]
        for sample_id in barcode_info:
            output_oriented_demultiplex_fastq = os.path.join(output_dir,oriented_basename) + "." + sample_id + ".fastq"
            if os.path.exists(output_oriented_demultiplex_fastq):
                output_sample_dir = os.path.join(output_dir,sample_id)
                create_dir(output_sample_dir)
                new_fastq_sample_path = os.path.join(output_sample_dir,oriented_basename) + "." + sample_id + ".fastq"
                shutil.move(output_oriented_demultiplex_fastq, new_fastq_sample_path)
                input_paths.append(new_fastq_sample_path)
                primerF.append(marker_info_primerF)
                primerR.append(marker_info_primerR)
                new_output_dirs.append(output_sample_dir)
                #Populate the dict_markers_outputDIR
                if "markerID" not in dict_markers_outputDIR.keys():
                    dict_markers_outputDIR["markerID"] = []
                    dict_markers_outputDIR["markerID"].append(output_sample_dir)
                elif "markerID" in dict_markers_outputDIR.keys():
                    dict_markers_outputDIR["markerID"].append(output_sample_dir)
        return input_paths,primerF,primerR,new_output_dirs,dict_markers_outputDIR
    else:
        log.info("An error has ocurred. Please check the oligos file - {}".format(oligos))
        return

#This function belongs to freq_table_uniqueness
#Delete all gaps (-) of a sequence which must be a string.
def remove_gaps(sequence):

    #Delete all gaps
    clean_seq = sequence.replace("-", "")
    return clean_seq

#Replace N for consensus_ambiguous_char selected by the user.
def use_set_ambiguous_char(sequence,consensus_ambiguous_char):

    #Convert to uppercase
    new_seq = sequence.upper()
    #Delete all gaps
    clean_seq = new_seq.replace("N",consensus_ambiguous_char)
    return clean_seq

#Process fastq to get the consensus
def process_fastq(input_fn, primer_a, primer_b, output_dir, primer_mismatch, min_base_score, min_seq_score, min_seq_cov,
                  max_len, aligner, sequence_max_mask, alignment_max_amb, max_len_delta,
                  expected_length, consensus_threshold, consensus_require_multiple,
                  mask_char, consensus_ambiguous_char, consensus_ignore_mask_char,min_SNP_cov,indel):
    # parse filename
    basename = os.path.splitext(os.path.basename(input_fn))[0]

    primer_b = str(Seq(primer_b).reverse_complement())
    
    # parse fastq file
    if max_len is None:
        records = (r for r in SeqIO.parse(input_fn, "fastq"))
    else:
        records = (r for r in SeqIO.parse(input_fn, "fastq") if len(r) < max_len)

    clean_records = list(trim_and_mask_seq_records(records, primer_a, primer_b, primer_mismatch,
                                                   min_base_score, basename, mask_char, min_seq_score))

    if len(clean_records) < min_seq_cov:
        log.info("alignment excluded - seq_count:{} < min_seq_cov:{} after trim_and_mask_seq_records - {}".format(
            len(clean_records), min_seq_cov, basename))
        return

    # define filter functions
    def mask_count_filter(r, _sequence_max_mask, _basename):
        _n_count = r.seq.upper().count(mask_char)
        if _n_count > _sequence_max_mask:
            log.info("seq excluded - mask_count:{} > sequence_max_mask:{} - {} {}".format(_n_count, _sequence_max_mask,
                                                                                          _basename, r.id))
            return False
        else:
            return True

    def len_variance_filter(r, _typical_len, _max_len_delta, _basename):
        len_delta = abs(len(r.seq) - _typical_len)
        if len_delta > _max_len_delta:
            log.info("seq excluded - len_delta:{} > max_len_delta:{} - {} {}".format(len_delta, _max_len_delta,
                                                                                     _basename, r.id))
            return False
        else:
            return True

    # apply mask_count_filter (based on coverage)
    if sequence_max_mask is not None:
        clean_records = [r for r in clean_records if mask_count_filter(r, sequence_max_mask, basename)]

    if len(clean_records) < min_seq_cov:
        log.info(
            "alignment excluded - seq_count:{} < min_seq_cov:{} after mask_count_filter - {}".format(
                len(clean_records),
                min_seq_cov, basename))
        return

    # determine typical_len for len_variance_filter
    # should be last per seq filter b/c requires stats on full seq list
    if max_len_delta is not None:
        if expected_length is not None:
            typical_len = expected_length
        else:
            try:
                typical_len = statistics.mode([len(r) for r in clean_records])
            except statistics.StatisticsError as _e:
                log.info("failover from mode to mean - {} - {}".format(_e, basename))
                try:
                    typical_len = statistics.mean([len(r) for r in clean_records])
                except statistics.StatisticsError as __e:
                    log.info("alignment excluded by length filter error - {} {} - {}".format(_e, __e, basename))
                    return
        # apply len variance filter
        clean_records = [r for r in clean_records
                         if len_variance_filter(r, typical_len, max_len_delta, basename)]

    if len(clean_records) < min_seq_cov:
        log.info(
            "alignment excluded - seq_count:{} < min_seq_cov:{} after len_variance_filter - {}".format(
                len(clean_records),
                min_seq_cov, basename))
        return


    # write trimmed and clean fasta files
    fasta_trim_cleaned_path = os.path.join(output_dir, basename) + "_Q" + str(min_base_score) + ".fasta"
    with open(fasta_trim_cleaned_path, "wt") as f:
        if len(clean_records) == 1:
            r = clean_records[0]
            mask_count = r.seq.upper().count(mask_char)
            r.description = "seq_count:1 mask_count:{} seq:len:{}".format(mask_count, len(r.seq))
            pass
        else:
            for r in clean_records:
                f.write(r.format("fasta"))

    # align trim and clean fasta
    number_seqs = len(clean_records)
    fasta_path = fasta_trim_cleaned_path
    output_path = os.path.splitext(fasta_path)[0] + ".aln.fasta"
    alignment_step(fasta_path,output_path,aligner,number_seqs,basename)
    input_consensus = output_path

    #Mask based on min coverage (minimum count)
    if min_SNP_cov:
        if len(clean_records) >= min_SNP_cov:
            aln_min_count_masked = mask_min_count_wraper(output_path,min_SNP_cov,mask_char,indel)
            #Convert aln to list of records
            if sequence_max_mask is not None:
                clean_records = [r for r in aln_min_count_masked if mask_count_filter(r, sequence_max_mask, basename)]
            else:
                clean_records = [r for r in aln_min_count_masked]

            if len(clean_records) < min_seq_cov:
                log.info(
                    "alignment excluded - seq_count:{} <= min_seq_cov:{} after mask_min_count_wraper - {}".format(
                    len(clean_records),
                    min_seq_cov, basename))
                return

            #Save the masked aligned
            fasta_trim_cleaned_min_cov_path = os.path.join(output_dir, basename) + "_Q" + str(min_base_score) + "_Cov" + str(min_SNP_cov) +".aln.fasta"
            with open(fasta_trim_cleaned_min_cov_path, "wt") as f:
                if len(clean_records) == 1:
                    r = clean_records[0]
                    mask_count = r.seq.upper().count(mask_char)
                    r.description = "seq_count:1 mask_count:{} seq:len:{}".format(mask_count, len(r.seq))
                    pass
                else:
                    for r in clean_records:
                        f.write(r.format("fasta"))
            input_consensus = fasta_trim_cleaned_min_cov_path    

        else:
            log.info(
                "alignment excluded - min_SNP_cov:{} <= number_seqs:{} after mask_min_cov - {}".format(
                    min_SNP_cov,len(clean_records),basename))
            return

    #Get a list with the path and frequency of each unique sequence. [[path1,number_seq1],[path2,number_seq2],...]
    new_fasta_path_freq_list = wraper_get_unique_seqs(input_consensus,basename,aligner,min_seq_cov,indel)

    #Iterate the unique sequences
    #Aligning unique sequences (haplotypes)
    for unique_sequences_fasta_path_freq in new_fasta_path_freq_list:
        number_seqs = unique_sequences_fasta_path_freq[1]
        basename = os.path.splitext(os.path.basename(unique_sequences_fasta_path_freq[0]))[0]
        if number_seqs != 1:
            if number_seqs >= min_seq_cov:
                #Prepare the paths and basename
                fasta_path = unique_sequences_fasta_path_freq[0]
                output_path = os.path.splitext(fasta_path)[0] + ".aln.fasta"

                #Realign unique sequences
                alignment_step(fasta_path,output_path,aligner,number_seqs,basename)

                input_consensus = output_path
                # write consensus fasta
                consensus = gap_consensus(input_consensus, threshold=consensus_threshold, mask_char=mask_char,
                              consensus_ambiguous_character="n",require_multiple=consensus_require_multiple,
                              consensus_ignore_mask_char=consensus_ignore_mask_char)
                amb_count = consensus.upper().count("N")
    
                #Check for the maximum number of ambiguous sites.
                if (alignment_max_amb is not None) and (amb_count > alignment_max_amb):
                    log.info("alignment of the haplotype excluded - amb_count:{} > alignment_max_amb:{} - {}".format(amb_count, alignment_max_amb,
                                                                                   basename))
                    continue
            
                consensus_seq = use_set_ambiguous_char(str(consensus),consensus_ambiguous_char)
                consensus_without_gaps = remove_gaps(consensus_seq)

                seq = Seq(data=consensus_without_gaps)
                description = "seq_count:{} amb_count:{} seq:len:{}".format(number_seqs, amb_count, len(consensus_without_gaps))
                # noinspection PyTypeChecker
                #Get the new name of the haplotype
                code_haplotype = basename.split(".")[-2].split("_")[-1]
                sampleID = os.path.basename(os.path.normpath(output_dir))
                haplotype_ID = sampleID + "_" + code_haplotype 
                seq_rec = SeqRecord(seq=seq, id=haplotype_ID , description=description)
                with open(os.path.join(output_dir, basename) + ".consensus.fasta", "wt") as f:
                    #fasta_entry = seq_rec.format("fasta").strip().split("\n")
                    #fasta_entry = fasta_entry[0] + "\n" + "".join(fasta_entry[1:]) + "\n"
                    f.write(seq_rec.format("fasta"))
            else:
                log.info("alignment of the haplotype excluded - min_SNP_cov:{} <= number_seqs:{} after resolve the haplotypes - {}".format(
                    min_SNP_cov,number_seqs,basename))
        elif number_seqs == 1:
            log.info("alignment of the haplotype excluded - number_seqs:{} after resolve the haplotypes  - {}".format(
                    number_seqs,basename))

    return


#Creates the final file with one or more of consensus sequence and the alignment. It also cleans up the intermediate files.  
def post_process(new_output_dirs,aligner,clean_up,min_prop_reads):
    list_final_consensus_files = []
    for output_dir in new_output_dirs:
        output_consensus_files = glob.glob(os.path.join(output_dir,'*consensus*'))
        sampleID = os.path.basename(os.path.normpath(output_dir))
        number_unique_sequences = len(output_consensus_files)
        
        if number_unique_sequences == 1:
            suffix_one_hap = "." + str(number_unique_sequences) + ".consensus.fasta"
            filename = sampleID + suffix_one_hap
            consensus_groups_fasta_path_one = os.path.join(output_dir,filename)
            list_final_consensus_files.append(consensus_groups_fasta_path_one)
            shutil.move(output_consensus_files[0],consensus_groups_fasta_path_one)

        elif number_unique_sequences > 1:
            #Divide the results in high frequency haplotypes (final_results_consensus) and minor_freq_haplotypes_consensus based on an arbitrary proportion of reads (min_prop_reads)
            final_results_consensus = []
            minor_freq_haplotypes_consensus = []

            dict_path_freq = {}
            number_reads_sum = 0
            
            #Get the frequency of each haplotype 
            for unique_consensus in output_consensus_files:
                number_reads = int(unique_consensus.split(".")[-3])
                dict_path_freq[unique_consensus] = number_reads
                number_reads_sum = number_reads_sum  + number_reads

            #Separate high and low frequency haplotypes
            for dict_path_freq_element in dict_path_freq.keys():
                porportion_reads = float(dict_path_freq[dict_path_freq_element])/float(number_reads_sum)
                if porportion_reads >= min_prop_reads:
                    final_results_consensus.append(dict_path_freq_element)
                else:
                    minor_freq_haplotypes_consensus.append(dict_path_freq_element)
            
            #Get the number of consensus per class
            len_final_results_consensus = len(final_results_consensus)
            len_minor_freq_haplotypes_consensus = len(minor_freq_haplotypes_consensus)
            
            if len_final_results_consensus > 0:
                #Get the file name
                suffix_freq = "." + str(len_final_results_consensus) + ".consensus.fasta"
                filename = sampleID + suffix_freq
                consensus_groups_fasta_path = os.path.join(output_dir,filename)
                
                #Save the results
                with open(consensus_groups_fasta_path, "w") as f:
                    for unique_consensus in final_results_consensus:
                        for record in SeqIO.parse(unique_consensus, "fasta"):
                            f.write(record.format("fasta"))
                list_final_consensus_files.append(consensus_groups_fasta_path)
                #Generate the alignment of consensus of unique sequences
                if len_final_results_consensus > 1:
                    output_path = os.path.splitext(consensus_groups_fasta_path)[0] + ".aln.fasta"
                    alignment_step(consensus_groups_fasta_path,output_path,aligner,number_unique_sequences,sampleID)

            if len_minor_freq_haplotypes_consensus > 0:
                #Get the file name
                suffix_low_freq = "." + str(len_minor_freq_haplotypes_consensus) + ".minor_freq_hap.consensus.fasta"
                filename = sampleID + suffix_low_freq
                consensus_minor_freq_groups_fasta_path = os.path.join(output_dir,filename)
                #Save the results
                with open(consensus_minor_freq_groups_fasta_path, "w") as f:
                    for unique_consensus in minor_freq_haplotypes_consensus:
                        for record in SeqIO.parse(unique_consensus, "fasta"):
                            f.write(record.format("fasta"))            

        if clean_up == "zip":
            #Compress the intermediate files
            all_files = glob.glob(os.path.join(output_dir,'*'))
            consensus_files = glob.glob(os.path.join(output_dir,'*consensus*'))
            
            if len(all_files) > 0:
                for consensus_file in consensus_files:
                    if consensus_file in all_files:
                        all_files.remove(consensus_file)

                intermediate_files_zip = sampleID + "_intermediate_files.zip"
                intermediate_files_zip_path = os.path.join(output_dir,intermediate_files_zip)
                with ZipFile(intermediate_files_zip_path, 'w') as zipObj2:
                    for file in all_files:
                        zipObj2.write(file,os.path.basename(file))
                        try:
                            os.remove(file)
                        except:
                            log.info("Error while deleting a file in post_process function: {}".format(file))

        elif clean_up == "remove":
            #Compress the intermediate files
            all_files = glob.glob(os.path.join(output_dir,'*'))
            consensus_files = glob.glob(os.path.join(output_dir,'*consensus*'))
            
            if len(all_files) > 0:
                for consensus_file in consensus_files:
                    if consensus_file in all_files:
                        all_files.remove(consensus_file)
            
                for file in all_files:
                    try:
                        os.remove(file)
                    except:
                        log.info("Error while deleting a file in post_process function: {}".format(file))
    return list_final_consensus_files

#Save all results of mode all in one file per marker
def consolidate_results(dict_markers_outputDIR,list_final_consensus_files,output_dir):
    final_results_dir = os.path.join(output_dir,"final_consensus_sequences")
    final_results_dir_marker_consensus_file_list = []
    create_dir(final_results_dir)
    if len(dict_markers_outputDIR.keys()) == 1 and list(dict_markers_outputDIR.keys())[0] == "markerID":
        final_results_file = os.path.join(final_results_dir,"consensus_sequences.fasta")
        final_results_dir_marker_consensus_file_list = [final_results_file]
        #Save the results into a file
        with open(final_results_file,"w") as final_results_file_handle:
            for file_result in list_final_consensus_files:
                with open(file_result,"r") as final_consensus_file_handle:
                    final_results_file_handle.write(final_consensus_file_handle.read())
        return final_results_dir_marker_consensus_file_list

    else:
        for markerID in dict_markers_outputDIR.keys():
            list_files_per_marker = []
            final_results_per_marker = os.path.join(final_results_dir,markerID)
            create_dir(final_results_per_marker)
            final_results_file = os.path.join(final_results_per_marker,markerID + ".fasta")
            final_results_dir_marker_consensus_file_list.append(final_results_file)
            #Save the consensus sequences in a file named as markerID.fasta
            with open(final_results_file,"w") as final_results_file_handle:
                for directory_results in dict_markers_outputDIR[markerID]:
                    for file_result in list_final_consensus_files:
                        if directory_results == os.path.dirname(file_result):
                            with open(file_result,"r") as final_consensus_file_handle:
                                final_results_file_handle.write(final_consensus_file_handle.read())
                                break
        return final_results_dir_marker_consensus_file_list

#check_file_exist belongs to create_results_info_file function
#Check the file
def check_file_exist(input_file):
    if os.path.isfile(input_file):
        if os.path.getsize(input_file) > 0:
            return "TRUE"
        else:
            return "FALSE"
    else:
        return "FALSE"

#create_table_oligos_file belongs to create_results_info_file function
#Process oligos file to get all sampleIDs per marker
def create_table_oligos_file(oligos):
    final_results_dict = {}
    markers = []
    samplesID = []
    with open(oligos, 'r') as file:
        reader = csv.reader(file,delimiter="\t")
        for line in reader:
            if line[0] == "primer" and len(line) == 3:
                markers.append("locusID")
            elif line[0] == "primer" and len(line) == 4:
                markers.append(line[3])
            elif line[0] == "barcode" and len(line) == 4:
                samplesID.append(line[3])
    for combination in itertools.product(samplesID,markers):
        #If only one marker, no markerID provided
        if combination[1] == "locusID":
            final_results_dict[combination[0]] = ["locusID",combination[0]]
        else:
            final_results_dict[combination[0] + "." + combination[1]] = [combination[1],combination[0]]
    return markers,final_results_dict

#process_mothur_demultiplex belongs to create_results_info_file function
#Process mothur demultiplex file to get the coverage per sample
def process_mothur_demultiplex(mothur_demultiplex,final_results_dict):
    demultiplexed_reads = []
    with open(mothur_demultiplex, 'r') as file:
        for line in file:
            sampleID = line.rstrip().split('\t')[1]
            demultiplexed_reads.append(sampleID)

    count_dict = Counter(demultiplexed_reads)
    for sampleID in final_results_dict.keys():
        if sampleID in count_dict.keys():
            final_results_dict[sampleID].append(str(count_dict[sampleID]))
        else:
            final_results_dict[sampleID].append("NA")
    return final_results_dict

#process_final_results belongs to create_results_info_file function
#Process final_consensus_sequences file to get information per haplotype
def process_final_results(final_results_dir_marker_consensus_file_list,final_results_dict,markers):
    final_results_list = []
    hap_info = []
    get_number_haplotype = []
    
    #Check when there is no markerID
    if markers == ["locusID"]:
        for final_consensus_unique_file in final_results_dir_marker_consensus_file_list:
            for record in SeqIO.parse(final_consensus_unique_file, "fasta"):
                haplotypeID = record.id
                sampleID = ("_").join(haplotypeID.split("_")[:-1])
                sampleID_info = [sampleID]
                sampleID_info = sampleID_info + final_results_dict[sampleID]
                sampleID_info.append(haplotypeID)
                seq_count = record.description.split(" ")[1].replace("seq_count:", "")
                sampleID_info.append(seq_count)
                num_ambiguity = record.description.split(" ")[2].replace("amb_count:", "")
                sampleID_info.append(num_ambiguity)
                seq_count = record.description.split(" ")[3].replace("seq:len:", "")
                sampleID_info.append(seq_count)
                #Save final result
                hap_info.append(sampleID_info)
                #Temporary list to get number of haplotypes
                get_number_haplotype.append(sampleID)
    
    else:
        for final_consensus_unique_file in final_results_dir_marker_consensus_file_list:
            #Get marker ID
            markerID = os.path.split(os.path.dirname(final_consensus_unique_file))[1]

            for record in SeqIO.parse(final_consensus_unique_file, "fasta"):
                haplotypeID = record.id
                sampleID = ("_").join(haplotypeID.split("_")[:-1]) + "." + markerID
                sampleID_info = [sampleID]
                sampleID_info = sampleID_info + final_results_dict[sampleID]
                sampleID_info.append(haplotypeID)
                seq_count = record.description.split(" ")[1].replace("seq_count:", "")
                sampleID_info.append(seq_count)
                num_ambiguity = record.description.split(" ")[2].replace("amb_count:", "")
                sampleID_info.append(num_ambiguity)
                seq_count = record.description.split(" ")[3].replace("seq:len:", "")
                sampleID_info.append(seq_count)
                #Save final result
                hap_info.append(sampleID_info)
                #Temporary list to get number of haplotypes
                get_number_haplotype.append(sampleID)
    
    #Obtain number of haplotypes per sample
    get_number_haplotype_dict = Counter(get_number_haplotype)
    get_number_haplotype_set = set(get_number_haplotype)
                
    #Add number of haplotypes per sample
    for idx in range(len(hap_info)):
        hap_info_line = hap_info[idx]
        number_haplotypes = str(get_number_haplotype_dict[hap_info_line[0]])
        hap_info[idx].insert(4,number_haplotypes)

    #Add missing data
    for sampleID in final_results_dict.keys():
        if sampleID not in get_number_haplotype_set:
            sample2add = [sampleID] + final_results_dict[sampleID] + ["NA"]*5
            hap_info.append(sample2add)

    return hap_info

#save_results_info_file belongs to create_results_info_file function
#Save the results into a file including header
def save_results_info_file(hap_info_info_final,results_info_file):
    with open(results_info_file,"w") as output_log_path_handle:
        #Save a header
        header = ["MarkerID","SampleID","Demultiplex_coverage","Number_haplotypes","HaplotypeID","Haplotype_coverage","Missing","Haplotype_length"]
        output_log_path_handle.write(("\t").join(header) + "\n")
        #Save results
        for hap_info in hap_info_info_final:
            new_hap_info = hap_info[1:]
            output_log_path_handle.write(("\t").join(new_hap_info) + "\n")
    return

#Process the results to get an info file
def create_results_info_file(output_oriented_fastq,final_results_dir_marker_consensus_file_list,oligos,results_info_file):
    #Get mothur demultiplex file
    mothur_demultiplex = os.path.splitext(output_oriented_fastq)[0] + ".group"
    #Check the files 
    if check_file_exist(mothur_demultiplex) == "FALSE":
        return
   
    for final_results_dir_marker_consensus_file in final_results_dir_marker_consensus_file_list:
        if check_file_exist(final_results_dir_marker_consensus_file) == "FALSE":
            return

    #Process oligos file to get all possible combinations markers and samples
    markers,final_results_dict = create_table_oligos_file(oligos)

    #Get coverage of demultiplexed samples
    final_results_dict = process_mothur_demultiplex(mothur_demultiplex,final_results_dict)

    #Get table with the information per sample
    hap_info_info_final = process_final_results(final_results_dir_marker_consensus_file_list,final_results_dict,markers)

    #Save into an output file
    save_results_info_file(hap_info_info_final,results_info_file)

    return

#save_results_info_file belongs to create_results_info_file_demultiplex function
#Save the results into a file including header
def save_results_info_file_demultiplex(final_results_dict,results_info_file):
    with open(results_info_file,"w") as output_log_path_handle:
        #Save a header
        header = ["MarkerID","SampleID","Demultiplex_coverage"]
        output_log_path_handle.write(("\t").join(header) + "\n")
        #Save results
        for sampleID in final_results_dict:
            sample_list = final_results_dict[sampleID]
            marker = sample_list[0].replace(sample_list[1] + ".","")
            new_sample_list = sample_list[1:]
            output_log_path_handle.write(marker + "\t" + ("\t").join(new_sample_list) + "\n")
    return

#Process the results to get an info file demultiplex
def create_results_info_file_demultiplex(output_oriented_fastq,oligos,results_info_file):
    #Get mothur demultiplex file
    mothur_demultiplex = os.path.splitext(output_oriented_fastq)[0] + ".group"
    #Check the files 
    if check_file_exist(mothur_demultiplex) == "FALSE":
        return

    #Process oligos file to get all possible combinations markers and samples
    markers,final_results_dict = create_table_oligos_file(oligos)

    #Get coverage of demultiplexed samples
    final_results_dict = process_mothur_demultiplex(mothur_demultiplex,final_results_dict)
    #Save into an output file
    save_results_info_file_demultiplex(final_results_dict,results_info_file)

    return

#Main function to get the consensus
def process_file_list(in_file_list, primer_mismatch, min_base_score, min_seq_score, min_seq_cov,
                      max_len, aligner, sequence_max_mask, alignment_max_amb,max_len_delta,
                      expected_length, consensus_threshold, consensus_require_multiple,
                      mask_char, consensus_ambiguous_char, consensus_ignore_mask_char,
                      min_SNP_cov,indel,number_processors=mp.cpu_count()):

    #Create the lists for infile and primers
    output_dirs=[]
    input_files=[]
    primer_F=[]
    primer_R=[]

    #Process the in_file_list
    with open(in_file_list, "r") as f:
        for l in f.readlines():
            line_to_list = l.strip().split("\t")
            if len(line_to_list) == 5:
                primer_F.append(line_to_list[1])
                primer_R.append(line_to_list[2])
                input_file_full_path = convert_to_full_path(line_to_list[3])
                input_files.append(input_file_full_path)
                output_dirs_full_path = convert_to_full_path(line_to_list[4])
                output_dirs.append(output_dirs_full_path)
            else:
                log.info("Check the in_file_list - {}".format(in_file_list))
                return
    
    if not len(output_dirs) == len(input_files) == len(primer_F) == len(primer_R):
        log.info("Check the in_file_list - {}".format(in_file_list))
        return

    #Create the dirs
    for output_dir in output_dirs:
        create_dir(output_dir)

    # create pool
    with mp.Pool(min(len(input_files), number_processors)) as p:
        p.starmap(process_fastq, zip(input_files,primer_F, primer_R,output_dirs,
                                     itertools.repeat(primer_mismatch),
                                     itertools.repeat(min_base_score),
                                     itertools.repeat(min_seq_score),
                                     itertools.repeat(min_seq_cov),
                                     itertools.repeat(max_len),
                                     itertools.repeat(aligner),
                                     itertools.repeat(sequence_max_mask),
                                     itertools.repeat(alignment_max_amb),
                                     itertools.repeat(max_len_delta),
                                     itertools.repeat(expected_length),
                                     itertools.repeat(consensus_threshold),
                                     itertools.repeat(consensus_require_multiple),
                                     itertools.repeat(mask_char),
                                     itertools.repeat(consensus_ambiguous_char),
                                     itertools.repeat(consensus_ignore_mask_char),
                                     itertools.repeat(min_SNP_cov),
                                     itertools.repeat(indel)
                                     ))
    return output_dirs


def main(_args):

    #Check possible errors in the supplied arguments
    check_args(_args.in_file_list,_args.in_file,_args.oligos,_args.aligner,_args.indel,_args.mode,_args.clean_up,_args.min_prop_reads)
    # call function based on batch or single mode

    #Convert to full paths
    if _args.in_file:
        _args.in_file = convert_to_full_path(_args.in_file)
    if _args.in_file_list:
        _args.in_file_list = convert_to_full_path(_args.in_file_list)
    if _args.oligos:
        _args.oligos = convert_to_full_path(_args.oligos)
    if _args.out_dir:
        _args.out_dir = convert_to_full_path(_args.out_dir)
    
    #Create a path for the summary results:
    summary_results = os.path.join(_args.out_dir,"summary_results")
    #Default of fixed parameters
    consensus_threshold = 0.7
    consensus_require_multiple = False

    #Process data according to the mode
    if _args.mode == "demultiplex":
        if _args.in_file is not None and _args.oligos is not None:
            output_dir = _args.out_dir

            # if output_path does not exist, create it
            log.info("Demultiplex mode - {}".format(_args.in_file))
            create_dir(output_dir)
            output_oriented_fastq,barcode_info,primer_names = demultiplex_wraper(_args.in_file,_args.oligos,output_dir,_args.number_processors,_args.primer_mismatch,_args.barcode_mismatch)
            
            #Create a summary file of demultiplex
            results_info_file = summary_results
            create_results_info_file_demultiplex(output_oriented_fastq,_args.oligos,results_info_file)
            return
        
        else:
            parser.print_help()
            sys.exit(1)

    elif _args.mode == "consensus":
        if _args.in_file_list is None:
            if _args.in_file is not None:
                log.info("Consensus mode - {}".format(_args.in_file))
                
                # if output_path does not exist, create it
                output_dir = _args.out_dir
                create_dir(output_dir)
                
                log.info("Demultiplex mode - {}".format(_args.in_file))
                process_fastq(_args.in_file, _args.primer_forward, _args.primer_reverse, output_dir, _args.primer_mismatch, _args.min_base_score, _args.min_seq_score,
                      _args.min_seq_cov,_args.max_len, _args.aligner,_args.sequence_max_mask,_args.alignment_max_amb,
                      _args.max_len_delta,_args.expected_length,consensus_threshold,consensus_require_multiple,
                      _args.mask_char, _args.consensus_ambiguous_char,_args.consensus_ignore_mask_char,_args.min_SNP_cov,_args.indel)
            
                clean_up = "keep"
                post_process([output_dir],_args.aligner,clean_up,_args.min_prop_reads)

                return

            else:
                parser.print_help()
                sys.exit(1)

    elif _args.mode == "all":
        if _args.in_file is not None and _args.oligos is not None:
            log.info("Demultiplex and consensus (all) mode - {}".format(_args.in_file))
        
            # if output_path does not exist, create it
            output_dir = _args.out_dir
            create_dir(output_dir)
        
            log.info("Demultiplex starts - {}".format(_args.in_file))        
            output_oriented_fastq,barcode_info,primer_names = demultiplex_wraper(_args.in_file,_args.oligos,output_dir,_args.number_processors,_args.primer_mismatch,_args.barcode_mismatch)

            if len(barcode_info) != 0:
                log.info("Demultiplex has successfully finished - {}".format(_args.in_file))
        
            input_paths,primerF,primerR,new_output_dirs,dict_markers_outputDIR = prepare_directories_files_all_mode(_args.in_file,output_dir,barcode_info,primer_names,_args.oligos)
        
            #Save a file for consensus_batch mode
        
            log.info("Consensus starts - {}".format(_args.in_file))
            consensus_batch_input = os.path.join(output_dir,"consensus_batch_input")
            log.info("Save an input file for consensus_batch mode - {}".format(consensus_batch_input))
            with open(consensus_batch_input,"w") as cbi_handler:
                for index in range(len(input_paths)):
                    cbi_handler.write("primer" + "\t" + primerF[index] + "\t" + primerR[index] + "\t" + input_paths[index] + "\t" + new_output_dirs[index] + "\n")
         
            #Run consensus 
            output_dirs = process_file_list(consensus_batch_input, _args.primer_mismatch, _args.min_base_score, _args.min_seq_score,
                          _args.min_seq_cov,_args.max_len,_args.aligner,_args.sequence_max_mask,_args.alignment_max_amb,
                          _args.max_len_delta,_args.expected_length,consensus_threshold,consensus_require_multiple,
                          _args.mask_char,_args.consensus_ambiguous_char,_args.consensus_ignore_mask_char,_args.min_SNP_cov,_args.indel,_args.number_processors)

            #Process the results
            list_final_consensus_files = post_process(output_dirs,_args.aligner,_args.clean_up,_args.min_prop_reads)
            final_results_dir_marker_consensus_file_list = consolidate_results(dict_markers_outputDIR,list_final_consensus_files,output_dir)
            
            results_info_file = summary_results
            create_results_info_file(output_oriented_fastq,final_results_dir_marker_consensus_file_list,_args.oligos,results_info_file)
            return

        else:
            parser.print_help()
            sys.exit(1)

    elif _args.mode == "consensus_batch":
        if _args.in_file_list is not None:
            log.info("Batch mode - {}".format(_args.in_file_list))
            #Run consensus
            output_dirs = process_file_list(_args.in_file_list, _args.primer_mismatch, _args.min_base_score, _args.min_seq_score,
                          _args.min_seq_cov,_args.max_len,_args.aligner,_args.sequence_max_mask,_args.alignment_max_amb,
                          _args.max_len_delta,_args.expected_length,consensus_threshold,consensus_require_multiple,
                          _args.mask_char,_args.consensus_ambiguous_char,_args.consensus_ignore_mask_char,_args.min_SNP_cov,_args.indel,_args.number_processors)

            clean_up = "keep"
            #Process the results
            list_final_consensus_files = post_process(output_dirs,_args.aligner,clean_up,_args.min_prop_reads)
            return

        else:
            parser.print_help()
            sys.exit(1)


    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":

    # argparse configuration
    # noinspection PyMissingOrEmptyDocstring
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser()

    # settings/options
    parser.add_argument('-l', '--aligner', help='The alignment software to use. Possible programs: clustalw, muscle or mafft. Default aligner: mafft', default="mafft")
    parser.add_argument('--mask_char', help='Character used to mask bases with quality below the threshold. Default: N', default="N")

    # files/directories
    parser.add_argument('-i', '--in_file', help='Path to input fastq file. This option will be used on demultiplex, consensus and all modes.')
    parser.add_argument('-o', '--out_dir', help='Path to output dir. This option will be used on demultiplex, consensus and all modes. Default output_dir: output.', default="output")
    #parser.add_argument('-out', '--out_results', help='Prefix of a tab delimited file summarizing the results. This option will be used on all mode. Default out_results: results.', default="results")
    parser.add_argument('-oligos', help='Path to oligos file for demultiplexing step (demultiplex and all modes). This is a tab delimited file with primers and barcodes following the pattern of oligos file of mothur program. For more information, please see https://mothur.org/wiki/oligos_file/.')

    #Cleanup
    parser.add_argument("--clean_up",help='Remove, compress or keep the intermediate files. This option is only activated with the all mode. Possible options: remove, zip or keep. Remove: eliminate all intermidiate files. zip: save the intermidiate files into a zip file. keep: save the intermediate files. Default: zip.',default="zip")

    #Mode
    parser.add_argument('--mode', help='Possible values: demultiplex, consensus, consensus_batch or all (demultiplex and consensus). Default: all.', default='all')
    # Consensus_batch mode filelist settings
    parser.add_argument('--in_file_list', help='Path to an input file for consensus_batch mode. This file must contain fields (no spaces allowed) separated by tab in the following order: Primer name, sequence of primer forward, sequence of primer reverse, path of the fastq file and path of output directory. Each line correponds to the information of one sample.', default=None)
    
    # Consensus mode primer sequences settings
    parser.add_argument('-F', '--primer_forward', help='Sequence of the primer forward. This option is required on the consensus mode.')
    parser.add_argument('-R', '--primer_reverse', help='Sequence of the primer reverse. This option is required on the consensus mode.')

    parser.add_argument('--number_processors', type=int, help='Set the number of processors. Mothur will be used all available processors for demultiplexing anyway. Total number of processors will be used as default.', default=mp.cpu_count())

    # base filters
    parser.add_argument('-n', '--min_base_score', type=int, help='Phred score below which a nucleotide is masked. Default: 30', default="30")

    # sequence filters
    parser.add_argument('-pm', '--primer_mismatch', type=int, help='Number of errors allowed in primer match for demultiplex, consensus, consensus_batch and all modes. This is equivalent to the flag pdiffs on mothur. Default: 2.', default="2")
    parser.add_argument('-bm', '--barcode_mismatch', type=int, help='Number of errors allowed in barcodes match for demultiplex mode. This is equivalent to the flag bdiffs on mother. Default: 2.', default="2")
    parser.add_argument('-p', '--sequence_max_mask', type=int, help='Number of mask characters allowed in the sequences. Sequences with higher number of bases masked will be removed. Default: 10.', default="10")
    parser.add_argument('-d', '--max_len_delta', type=int, help='Maximum allowed discrepancy from the length mode of the sequences. Default: 5.', default="5")
   
    # optional sequence filters
    parser.add_argument('--expected_length', type=int, help='Use a user specified length instead of the modal length, when applying the max_len_delta filter. This flag is optional.', default=None)
    parser.add_argument('--max_len', type=int, help='Optional, max length overwhich a sequence is excluded.', default=None)
    parser.add_argument('--min_seq_score', type=int, help='Optional, average phred score below which a sequence is excluded.', default=None)
    
    # Optional variant calling filters
    parser.add_argument('--min_SNP_cov', type=int, help='Minimum coverage (number of counts) to consider a single nucleotide polymorphism (SNP). In case indel option is T, indels will be trated as SNPs. All polimorphic bases with less than min_SNP_cov will be masked. Default: 5', default="5")
    parser.add_argument('--indel', help='Treatment of indels (-) to split up sequences into haplotypes. Options NS (new state), N (nucleotide), and M (missing data). NS: gaps are treated as a new state (a character different than a nucleotide). N: gaps are treated as a fifth state (equivalent to a nucleotide), M: gaps are treated as missing data. Default: M', default='M')

    # alignment filters
    parser.add_argument('-s', '--min_seq_cov', type=int, help='Minimum coverage (counts) of sequences before alignment is excluded. Default: 10.', default="10")
    parser.add_argument('-t', '--alignment_max_amb', type=int, help='Number of consensus_ambiguous_char allowed in final consensus sequences. Default: 0.', default="0")
    parser.add_argument('--min_prop_reads',type=float, help='Minimum proportion of reads to call an haplotype. Consensus sequences with less than this proportion (frequency) will be stored in another file with a suffix: .minor_freq_hap.consensus.fasta. Default: 0.2.', default="0.2")

    # consensus options
    #parser.add_argument('--consensus_threshold', type=float, help='Proportion threshold for consensus to call a base. Default: 0.7.', default=".7")
    #parser.add_argument('--consensus_require_multiple', help='Require multiple instanses for consensus to call a base.', action='store_true')
    parser.add_argument('--consensus_ambiguous_char', help='Character representing bases under consensus threshold. Default: N.', default="N")
    parser.add_argument('--consensus_ignore_mask_char', help='Discount mask character when calculating consensus.', action='store_true')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # run main
    try:
        exit(main(args))
    except Exception as e:
        log.exception("Exception in main(): {}".format(e))
        exit(1)
