#!/usr/bin/env python
# coding=utf-8
"""
If running as a script, run with the -h flag for usage documentation.
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
from Bio import pairwise2
from Bio.Align import AlignInfo
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio.Alphabet import IUPAC
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

log = logging.getLogger("ccs-consensuser.py")


# some functions from: https://github.com/chapmanb/bcbb/blob/master/align/adaptor_trim.py
# _remove_adaptor
# trim_adaptor
# trim_adaptor_w_qual

def check_args(in_file_list,in_file,oligos,aligner,indel,mode,clean_up):
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
            log.info("Simple spaces are not allowed, please use tab as field separator - {}".format(_args.in_file_list))
            sys.exit(1)
        elif result == "Path":
            log.info("The file does not exist - {}".format(_args.in_file_list))
            sys.exit(1)           
    if in_file is not None:
        result = check_whitespaces(in_file)
        if result == "whitespace":
            log.info("Simple spaces are not allowed, please use tab as field separator - {}".format(in_file))
            sys.exit(1)
        elif result == "Path":
            log.info("The file does not exist - {}".format(in_file))
            sys.exit(1)
    if oligos is not None:
        result = check_whitespaces(oligos)
        if result == "whitespace":
            log.info("Simple spaces are not allowed, please use tab as field separator - {}".format(oligos))
            sys.exit(1)
        elif result == "Path":
            log.info("The file does not exist - {}".format(oligos))
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

    if indel not in ["Y","N"]:
        log.info("indel must be Y or N")
        sys.exit(1)

    if clean_up not in ["remove","zip","keep"]:
        log.info("clean_up must be remove, zip or keep.")
        sys.exit(1)        

    return

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
                primer_names.append(["sampleID",line[1],line[2]])
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
        barcode_name_F = "primerF" + str(count)
        barcode_name_R = "primerR" + str(count)
        primer_info_new.append([barcode_name_F,primer[0]])
        primer_info_new.append([barcode_name_R,primer[1]])
        primerF.append(barcode_name_F)
        primerR.append(barcode_name_R)
        count += 1

    with open(primer_fasta_path, 'w') as primer_fasta:
        for primer in primer_info_new:
            primer_fasta.write(">" + primer[0] + "\n" + primer[1] + "\n")
    return primerF,primerR,barcode_info,primer_names,primer_fasta_path

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
    #makeblastdb -in its_anas.fasta -input_type fasta -dbtype nucl -out its_anas.fasta
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
        #blastn -task blastn-short -db its_anas.fasta -query primer  -out its_anas.fasta_primer -evalue 1e-1 -outfmt 6 -max_target_seqs 2500  -num_threads 1
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
        os.remove(blast_db + ".nhr")
        os.remove(blast_db + ".nin")
        os.remove(blast_db + ".nsq")
    
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
        #blastn -task blastn-short -db its_anas.fasta -query primer  -out its_anas.fasta_primer -evalue 1e-1 -outfmt 6 -max_target_seqs 2500  -num_threads 1
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

#This function re-orient the reads and runs the fastq.info of the mothur program to demultiplex the reads.
def demultiplex_wraper(input_fastq,oligos_path,output_dir,number_processors,primer_mismatch,barcode_mismatch):
    list_primerF,list_primerR,barcode_info,primer_names,primer_fasta_path = parse_oligos(oligos_path,output_dir)
    input_fasta,number_seq = convert_fastq2fasta(input_fastq,output_dir)
    output_blast = short_blast(input_fasta,output_dir,number_seq,primer_fasta_path,number_processors)
    orientation_dict = convert_blastoutfmt6_2_dict(output_blast,list_primerF,list_primerR)
    output_oriented_fastq = reorientation(input_fastq,output_dir,orientation_dict)
    run_mothur(output_oriented_fastq,oligos_path,primer_mismatch,barcode_mismatch)
    return barcode_info,primer_names

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
    exact_pos = str(seq).find(adaptor)
    if exact_pos >= 0:
        seq_region = str(seq[exact_pos:exact_pos + len(adaptor)])
        adapt_region = adaptor
    else:
        aligns = pairwise2.align.localms(str(seq), str(adaptor),
                                         5.0, -4.0, -9.0, -0.5, one_alignment_only=True,
                                         gap_char=gap_char)
        if len(aligns) == 0:
            adapt_region, seq_region = ("", "")
        else:
            seq_a, adaptor_a, score, start, end = aligns[0]
            adapt_region = adaptor_a[start:end]
            seq_region = seq_a[start:end]
    matches = match_func_amb(seq_region,adapt_region)
    # too many errors -- no trimming
    if (len(adaptor) - matches) > primer_mismatch:
        return seq
    # remove the adaptor sequence and return the result
    else:
        return _remove_adaptor(seq, seq_region.replace(gap_char, ""),
                               right_side)

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


def gap_consensus(input_consensus, threshold=.7, mask_char="N", consensus_ambiguous_char="X",
                  consensus_alpha=None, require_multiple=False, consensus_ignore_mask_char=False):
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
        alignment = AlignIO.read(f, "fasta", alphabet=IUPAC.ambiguous_dna)

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
            consensus += consensus_ambiguous_char
        elif (len(max_atoms) == 1) and ((float(max_size) / float(num_atoms)) >= threshold):
            consensus += max_atoms[0]
        else:
            consensus += consensus_ambiguous_char

    # we need to guess a consensus alphabet if one isn't specified
    if consensus_alpha is None:
        # noinspection PyProtectedMember
        consensus_alpha = summary_align._guess_consensus_alphabet(consensus_ambiguous_char)

    return Seq(consensus, consensus_alpha)


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

    masked_seq_req.seq = Seq("".join(base_list), alphabet=IUPAC.ambiguous_dna)

    return masked_seq_req


def trim_and_mask_seq_records(records, primer_a, primer_b, primer_mismatch, min_base_score, basename, mask_char,
                              min_seq_score=None):
    for seq_rec in records:
        trimed_seq_rec = trim_both_ends(seq_rec, primer_a, primer_b, primer_mismatch, reverse_complement=False)

        # primers found in forward direction
        if trimed_seq_rec is not None:
            avg_score = np.mean(trimed_seq_rec.letter_annotations["phred_quality"])
            if min_seq_score and (avg_score < min_seq_score):
                log.info("seq excluded - avg_score:{:4.2f} < min_seq_score:{} - {} {}".format(avg_score, min_seq_score,
                                                                                              basename, seq_rec.id))
                continue
            yield mask_seq_record(trimed_seq_rec, min_base_score, mask_char=mask_char, inplace=True)

        # primers not found in forward direction
        else:
            trimed_seq_rec = trim_both_ends(seq_rec, primer_a, primer_b, primer_mismatch, reverse_complement=True)
            if trimed_seq_rec is not None:  # primers found in reverse direction
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


def alignment_step(fasta_path,output_path,aligner,number_seqs,basename):
    
    if aligner == "muscle":
        cline = MuscleCommandline(input=fasta_path,
                                  out=output_path, fasta=True)
        try:
            # noinspection PyUnusedLocal
            stdout, stderr = cline()
            log.debug(stderr)
        except Bio.Application.ApplicationError as _e:
            log.info("alignment failed - {} - {}".format(_e, basename))
            return

    elif aligner == "clustalw":
        cline = ClustalwCommandline("clustalw2", infile=fasta_path,
                                  outfile=output_path, output="FASTA")
        try:
            # noinspection PyUnusedLocal
            stdout, stderr = cline()
            log.debug(stderr)
        except Bio.Application.ApplicationError as _e:
            log.info("alignment failed - {} - {}".format(_e, basename))
            return

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
        
        #logging.info(clusterID + ": mafft alignment has finished the analysis.")
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
def mask_min_count_wraper(aln_fasta,min_cov,mask_char,indel):
    
    aln_dna = AlignIO.read(open(aln_fasta), 'fasta')

    dict_pos={}
    for python_pos in range(aln_dna.get_alignment_length()):
        pos_aln = aln_dna[:,python_pos].upper()
        base_counter = Counter(pos_aln)
        low_coverage = analize_counter(base_counter,min_cov,indel)
        if len(low_coverage) > 0:
            dict_pos[python_pos] = low_coverage

    records = mask_low_coverage(aln_dna,mask_char,dict_pos)

    return records

#Return the characters with low coverage (based on min_cov). Belongs to mask_min_count_wraper(aln_fasta,min_cov).
def analize_counter(base_counter,min_cov,indel):
    low_coverage = []
    bases = ["A","T","C","G"]
    if indel == "Y":
        bases.append("-")
    for base in base_counter.keys():
        if base_counter[base] < min_cov and base in bases:
            low_coverage.append(base)

    return low_coverage

#Mask the bases based on a dictionary (key = python position and values list of bases to be masked). Belongs to mask_min_count_wraper(aln_fasta,min_cov).
def mask_low_coverage(aln_dna,mask_char,dict_pos):
    for idx_record in range(len(aln_dna)):
        seq = str(aln_dna[idx_record].seq)
        #for pos_low_coverage in dict_pos.keys():
        #   if seq[pos_low_coverage] in dict_pos[pos_low_coverage]:
        #       seq[pos_low_coverage] = "N"
        
        new_seq = seq[:]
        for base_pos in range(len(seq) - 1, -1, -1):
            if base_pos in dict_pos.keys() and seq[base_pos] in dict_pos[base_pos]:
                new_seq = new_seq[:base_pos] + mask_char + new_seq[base_pos+1:]
        #   elif seq[base_pos] == "-":
        #       new_seq = new_seq[:base_pos] + "" + new_seq[base_pos+1:]
        aln_dna[idx_record].seq = Seq(new_seq)

    return aln_dna

#Wraper to get the "haplotypes" (unique sequences). It returns a list of path to fasta unaligned sequences and haplotype frequency (number of seqs per file).
def wraper_get_unique_seqs(path_fasta,basename,aligner,min_seq_count,indel):
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

    new_fasta_path_freq_list = dict_unique_to_files(path_fasta, basename, dict_unique,min_seq_count)

    return new_fasta_path_freq_list

#This function belongs to wraper_get_unique_seqs
#Compare two sequences and return the number of matches between sequences. Missing data (N) is considered a wildcard, so it will match any base or even a gap.
def compare_seq(seq1,seq2,indel):
    counter = 0
    if indel == "N":
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
    elif indel == "Y":
        for base_idx in range(len(seq1)):
            if seq1[base_idx] == seq2[base_idx]:
                counter = counter + 1
                continue
            elif seq1[base_idx] == "N" or seq2[base_idx] == "N":
                counter = counter + 1
                continue
    return counter

#This function belongs to wraper_get_unique_seqs
#Convert a dictionary with haplotypes (id1:[id1,id2,...]) to separate files.
def dict_unique_to_files(path_fasta_aln, basename, dict_unique,min_seq_count):
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
    if primer_names[0][0] != "sampleID":
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

        #Remove empty marker folders
        for marker_info in primer_names:
            marker_directory = os.path.join(output_dir,marker_info[0])
            if len(os.listdir(marker_directory)) == 0:
                os.rmdir(marker_directory)    
        return input_paths,primerF,primerR,new_output_dirs
    
    #For one primer pair (markers)
    elif len(primer_names) == 1 and primer_names[0][0] == "sampleID":
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
        return input_paths,primerF,primerR,new_output_dirs
    else:
        log.info("An error has ocurred. Please check the oligos file - {}".format(oligos))
        return

#This function belongs to freq_table_uniqueness
#Delete all gaps (-) of a sequence which must be a string.
def remove_gaps(sequence):

    #Convert to uppercase
    new_seq = sequence.upper()
    #Delete all gaps
    clean_seq = new_seq.replace("-", "")
    return clean_seq


#To do:Add min_count and gap=True options
'''
Default:input_fn, primer_a, primer_b, output_dir, primer_mismatch, min_base_score, min_seq_score=None, min_seq_count=1,
                  max_len=None, aligner, sequence_max_mask=None, alignment_max_amb=None, max_len_delta=None,
                  expected_length=None, consensus_threshold=0.7, consensus_require_multiple=False,
                  mask_char="N", consensus_ambiguous_char="X", consensus_ignore_mask_char=False,min_cov,indel):
'''
def process_fastq(input_fn, primer_a, primer_b, output_dir, primer_mismatch, min_base_score, min_seq_score, min_seq_count,
                  max_len, aligner, sequence_max_mask, alignment_max_amb, max_len_delta,
                  expected_length, consensus_threshold, consensus_require_multiple,
                  mask_char, consensus_ambiguous_char, consensus_ignore_mask_char,min_cov,indel):
    # parse filename
    basename = os.path.splitext(os.path.basename(input_fn))[0]

    primer_b = str(Seq(primer_b).reverse_complement())
    
    # parse fastq file
    if max_len is None:
        records = (r for r in SeqIO.parse(input_fn, "fastq", alphabet=IUPAC.ambiguous_dna))
    else:
        records = (r for r in SeqIO.parse(input_fn, "fastq", alphabet=IUPAC.ambiguous_dna) if len(r) < max_len)

    clean_records = list(trim_and_mask_seq_records(records, primer_a, primer_b, primer_mismatch,
                                                   min_base_score, basename, mask_char, min_seq_score))

    if len(clean_records) < min_seq_count:
        log.info("alignment excluded - seq_count:{} < min_seq_count:{} after trim_and_mask_seq_records - {}".format(
            len(clean_records), min_seq_count, basename))
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

    if len(clean_records) < min_seq_count:
        log.info(
            "alignment excluded - seq_count:{} < min_seq_count:{} after mask_count_filter - {}".format(
                len(clean_records),
                min_seq_count, basename))
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
                    log.info("alignment excluded - {} {} - {}".format(_e, __e, basename))
                    return
        # apply len variance filter
        clean_records = [r for r in clean_records
                         if len_variance_filter(r, typical_len, max_len_delta, basename)]

    if len(clean_records) < min_seq_count:
        log.info(
            "alignment excluded - seq_count:{} < min_seq_count:{} after len_variance_filter - {}".format(
                len(clean_records),
                min_seq_count, basename))
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
    #TO DO list if dont pass due to low number of sequences save the file basename.number_seq.fastq
    if min_cov:
        if len(clean_records) >= min_cov:
            aln_min_count_masked = mask_min_count_wraper(output_path,min_cov,mask_char,indel)
            #Convert aln to list of records
            if sequence_max_mask is not None:
                clean_records = [r for r in aln_min_count_masked if mask_count_filter(r, sequence_max_mask, basename)]
            else:
                clean_records = [r for r in aln_min_count_masked]

            if len(clean_records) < min_seq_count:
                log.info(
                    "alignment excluded - seq_count:{} <= min_seq_count:{} after mask_min_count_wraper - {}".format(
                    len(clean_records),
                    min_seq_count, basename))
                return

            #Save the masked aligned
            fasta_trim_cleaned_min_cov_path = os.path.join(output_dir, basename) + "_Q" + str(min_base_score) + "_Cov" + str(min_cov) +".aln.fasta"
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
                "alignment excluded - min_cov:{} <= number_seqs:{} after mask_min_cov - {}".format(
                    min_cov,len(clean_records),basename))
            return

    #Get a list with the path and frequency of each unique sequence. [[path1,number_seq1],[path2,number_seq2],...]
    new_fasta_path_freq_list = wraper_get_unique_seqs(input_consensus,basename,aligner,min_seq_count,indel)

    #Iterate the unique sequences
    #Aligning unique sequences (haplotypes)
    for unique_sequences_fasta_path_freq in new_fasta_path_freq_list:
        number_seqs = unique_sequences_fasta_path_freq[1]

        if number_seqs != 1 and number_seqs >= min_seq_count:
            #Prepare the paths and basename
            fasta_path = unique_sequences_fasta_path_freq[0]
            basename = os.path.splitext(os.path.basename(fasta_path))[0]
            output_path = os.path.splitext(fasta_path)[0] + ".aln.fasta"

            #Realign unique sequences
            alignment_step(fasta_path,output_path,aligner,number_seqs,basename)

            input_consensus = output_path
            # write consensus fasta
            consensus = gap_consensus(input_consensus, threshold=consensus_threshold, mask_char=mask_char,
                              consensus_ambiguous_char=consensus_ambiguous_char, consensus_alpha=IUPAC.ambiguous_dna,
                              require_multiple=consensus_require_multiple,
                              consensus_ignore_mask_char=consensus_ignore_mask_char)
            amb_count = consensus.upper().count(mask_char)
    
            #Check for the maximum number of ambiguous sites.
            if (alignment_max_amb is not None) and (amb_count > alignment_max_amb):
                log.info(
                    "alignment excluded - amb_count:{} > alignment_max_amb:{} - {}".format(amb_count, alignment_max_amb,
                                                                                   basename))
                return
            seq = Seq(data=str(consensus), alphabet=IUPAC.ambiguous_dna)
            description = "seq_count:{} amb_count:{} seq:len:{}".format(number_seqs, amb_count, len(consensus))
            # noinspection PyTypeChecker
            seq_rec = SeqRecord(seq=seq, id=basename , description=description)
            with open(os.path.join(output_dir, basename) + ".consensus.fasta", "wt") as f:
                fasta_entry = seq_rec.format("fasta").strip().split("\n")
                fasta_entry = fasta_entry[0] + "\n" + "".join(fasta_entry[1:]) + "\n"
                f.write(fasta_entry)

def post_process(new_output_dirs,aligner,clean_up):
    for output_dir in new_output_dirs:
        output_consensus_files = glob.glob(os.path.join(output_dir,'*consensus*'))
        sampleID = os.path.basename(os.path.normpath(output_dir))
        number_unique_sequences = len(output_consensus_files)
        suffix = "." + str(number_unique_sequences) + ".consensus.fasta"
        filename = sampleID + suffix
        consensus_groups_fasta_path = os.path.join(output_dir,filename)
        
        if number_unique_sequences == 1:
            shutil.move(output_consensus_files[0],consensus_groups_fasta_path)

        elif number_unique_sequences > 1:
            with open(consensus_groups_fasta_path, "w") as f:
                for unique_consensus in output_consensus_files:
                    for record in SeqIO.parse(unique_consensus, "fasta"):
                        sequence = str(record.seq)
                        new_sequence = remove_gaps(sequence)
                        record.seq = Seq(new_sequence)
                        f.write(record.format("fasta"))
        
            #Generate the alignment of consensus of unique sequences
            output_path = os.path.splitext(consensus_groups_fasta_path)[0] + ".aln.fasta"
            alignment_step(consensus_groups_fasta_path,output_path,aligner,number_unique_sequences,sampleID)

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
                        zipObj2.write(file)
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
    return


'''
Default parameters:
in_file_list, primer_mismatch, min_base_score, min_seq_score=None, min_seq_count=1,max_len=None, 
aligner=mafft, sequence_max_mask=None, alignment_max_amb=None,max_len_delta=None,expected_length=None, 
consensus_threshold=0.7, consensus_require_multiple=False,mask_char="N", consensus_ambiguous_char="X", 
consensus_ignore_mask_char=False,min_cov,indel="Y",number_processors=mp.cpu_count()

'''

def process_file_list(in_file_list, primer_mismatch, min_base_score, min_seq_score, min_seq_count,
                      max_len, aligner, sequence_max_mask, alignment_max_amb,max_len_delta,
                      expected_length, consensus_threshold, consensus_require_multiple,
                      mask_char, consensus_ambiguous_char, consensus_ignore_mask_char,
                      min_cov,indel,number_processors=mp.cpu_count()):

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
                input_files.append(line_to_list[3])
                output_dirs.append(line_to_list[4])
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
                                     itertools.repeat(min_seq_count),
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
                                     itertools.repeat(min_cov),
                                     itertools.repeat(indel)
                                     ))
    return output_dirs


def main(_args):
    #Check possible errors in the supplied arguments
    check_args(_args.in_file_list,_args.in_file,_args.oligos,_args.aligner,_args.indel,_args.mode,_args.clean_up)
    # call function based on batch or single mode

    if _args.mode == "demultiplex":
        if _args.in_file is not None and _args.oligos is not None:
            # if output_path does not exist, create it
            log.info("Demultiplex mode - {}".format(_args.in_file))
            output_dir = _args.out_dir
            create_dir(output_dir)
            demultiplex_wraper(_args.in_file,_args.oligos,output_dir,_args.number_processors,_args.primer_mismatch,_args.barcode_mismatch)
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
                      _args.min_seq_count,_args.max_len, _args.aligner,_args.sequence_max_mask,_args.alignment_max_amb,
                      _args.max_len_delta,_args.expected_length,_args.consensus_threshold,_args.consensus_require_multiple,
                      _args.mask_char, _args.consensus_ambiguous_char,_args.consensus_ignore_mask_char,_args.min_cov,_args.indel)
            
                post_process([output_dir],_args.aligner,_args.clean_up)

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

            barcode_info,primer_names = demultiplex_wraper(_args.in_file,_args.oligos,output_dir,_args.number_processors,_args.primer_mismatch,_args.barcode_mismatch)

            if len(barcode_info) != 0:
                log.info("Demultiplex has successfully finished - {}".format(_args.in_file))
        
            input_paths,primerF,primerR,new_output_dirs = prepare_directories_files_all_mode(_args.in_file,output_dir,barcode_info,primer_names,_args.oligos)
        
            #Save a file for consensus_batch mode
        
            log.info("Consensus starts - {}".format(_args.in_file))
            consensus_batch_input = os.path.join(output_dir,"consensus_batch_input")
            log.info("Save an input file for consensus_batch mode - {}".format(consensus_batch_input))
            with open(consensus_batch_input,"w") as cbi_handler:
                for index in range(len(input_paths)):
                    cbi_handler.write("primer" + "\t" + primerF[index] + "\t" + primerR[index] + "\t" + input_paths[index] + "\t" + new_output_dirs[index] + "\n")
         
            #Run conensus 
            output_dirs = process_file_list(consensus_batch_input, _args.primer_mismatch, _args.min_base_score, _args.min_seq_score,
                          _args.min_seq_count,_args.max_len,_args.aligner,_args.sequence_max_mask,_args.alignment_max_amb,
                          _args.max_len_delta,_args.expected_length,_args.consensus_threshold,_args.consensus_require_multiple,
                          _args.mask_char,_args.consensus_ambiguous_char,_args.consensus_ignore_mask_char,_args.min_cov,_args.indel,_args.number_processors)       
            
            post_process(output_dirs,_args.aligner,_args.clean_up)

            return

        else:
            parser.print_help()
            sys.exit(1)

    elif _args.mode == "consensus_batch":
        if _args.in_file_list is not None:
            log.info("Batch mode - {}".format(_args.in_file_list))
        
            output_dirs = process_file_list(_args.in_file_list, _args.primer_mismatch, _args.min_base_score, _args.min_seq_score,
                          _args.min_seq_count,_args.max_len,_args.aligner,_args.sequence_max_mask,_args.alignment_max_amb,
                          _args.max_len_delta,_args.expected_length,_args.consensus_threshold,_args.consensus_require_multiple,
                          _args.mask_char,_args.consensus_ambiguous_char,_args.consensus_ignore_mask_char,_args.min_cov,_args.indel,_args.number_processors) 
            
            post_process(output_dirs,_args.aligner,_args.clean_up)
            return

        else:
            parser.print_help()
            sys.exit(1)


    else:
        parser.print_help()
        sys.exit(1)

#To do list
#Option to clean intermediate files
#align the consensus
#Deal with mothur.1635898673.logfile (input directory)

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
    parser.add_argument('--mask_char', help='Character used to mask bases with quality below threshold. Default: N', default="N")

    # files/directories
    parser.add_argument('-i', '--in_file', help='Path to input fastq file. This option will be used on demultiplex, consensus and all modes.')
    parser.add_argument('-o', '--out_dir', help='Path to output dir. This option will be used on demultiplex, consensus and all modes. Default output_dir: output.', default="output")
    parser.add_argument('-oligos', help='Path to oligos file for demultiplexing step (demultiplex and all modes). This is a tab delimited file with primers and barcodes following the pattern of oligos file of mothur program. For more information, please see https://mothur.org/wiki/oligos_file/.')
    parser.add_argument('-F', '--primer_forward', help='Sequence of the sequence of the primer forward. It will be used on consensus mode.')
    parser.add_argument('-R', '--primer_reverse', help='Sequence of the sequence of the primer reverse. It will be used on consensus mode.')
    #parser.add_argument('--overwrite', help='Set this flag to disable creating new out dir', action='store_true')

    #Cleanup
    parser.add_argument("--clean_up",help='Remove, compress or keep the intermediate files. Possible options: remove, zip or keep. Default: zip.',default="zip")

    #Mode
    parser.add_argument('--mode', help='Possible values: demultiplex, consensus, consensus_batch or all (demultiplex and consensus). Default: all.', default='all')
    # batch mode filelist settings
    parser.add_argument('--in_file_list', help='Path to text to be used on consensus_batch mode. Fields separated by tab are: Primer name, sequence of primer forward, sequence of primer reverse, path of the fastq file and path of output directory. The information of each sample must be on each line.', default=None)
    #parser.add_argument('--in_dir', help='input dir of files named in in_file_list', default=None)
    parser.add_argument('--number_processors', type=int, help='Set the number of processors. Mothur will be used all available processors for demultiplexing anyway. Total number of processors will be used as default.', default=mp.cpu_count())

    # base filters
    parser.add_argument('-n', '--min_base_score', type=int, help='Phred score below which a nucleotide is masked. Default: 30', default="30")

    # sequence filters
    parser.add_argument('-pm', '--primer_mismatch', type=int, help='Number of errors allowed in primer match. For demultiplex (mothur flag pdiffs), consensus, consensus_batch and all modes. Default: 2.', default="2")
    parser.add_argument('-bm', '--barcode_mismatch', type=int, help='Number of errors allowed in barcodes match. For demultiplex mode (mothur flag bdiffs). Default: 2.', default="2")
    parser.add_argument('-p', '--sequence_max_mask', type=int, help='Number of mask_char allowed in sequences to be aligned. Default: 5.', default="5")
    parser.add_argument('-d', '--max_len_delta', type=int, help='Allowed maximum variation from mode of sequence length. Default: 5.',
                        default="5")
    # Optional variant calling
    parser.add_argument('--min_cov', type=int, help='Minimum coverage (number of counts) to accept a variant (SNP or indel). Otherwise the variant will be masked. Default: 0', default="0")
    parser.add_argument('--indel', help='Consider an indel as a polymorphic site. This option must be used along with min_cov. Possible options Y or N. WARNING. Indels might be masked. Default: N', default='N')

    # optional sequence filters
    parser.add_argument('--expected_length', type=int, help='Optional, replaces average in max_len_delta filter.', default=None)
    parser.add_argument('--max_len', type=int, help='Optional, max length overwhich a sequence is excluded.', default=None)
    parser.add_argument('--min_seq_score', type=int, help='Optional, average phred score below which a sequence is excluded.', default=None)

    # alignment filters
    parser.add_argument('-s', '--min_seq_count', type=int, help='Minimum count of sequences before alignment is excluded. Default: 5.',
                        default="5")
    parser.add_argument('-t', '--alignment_max_amb', type=int, help='Number of consensus_ambiguous_char allowed in final consensus sequences. Default: 5.', default="5")

    # consensus options
    parser.add_argument('--consensus_threshold', type=float, help='Proportion threshold for consensus to call a base. Default: 0.7.',
                        default=".7")
    parser.add_argument('--consensus_require_multiple', help='Require multiple instanses for consensus to call a base.',
                        action='store_true')
    parser.add_argument('--consensus_ambiguous_char', help='Character representing bases under consensus threshold. Default: n.',
                        default="n")
    parser.add_argument('--consensus_ignore_mask_char', help='Discount mask character when calculating consensus.',
                        action='store_true')

    args = parser.parse_args()

    # run main
    try:
        exit(main(args))
    except Exception as e:
        log.exception("Exception in main(): {}".format(e))
        exit(1)
