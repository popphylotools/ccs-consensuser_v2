from Bio import SeqIO
from Bio.Seq import Seq
import os,sys,argparse
from glob import glob
import csv
import random
import copy
import subprocess
import itertools
import multiprocessing as mp


"""
FUNCTIONS
"""
#Convert a relative path to full path
def convert_to_full_path(path):
    absolute_path = os.path.abspath(path)
    return absolute_path

#Create a directory
def create_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    return dir_path

#Convert fasta file to a list of id and seq
def parse_fasta(input_fasta):
    list_records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        list_records.append(record)

    if len(list_records) == 1:
        id = list_records[0].id
        seq = str(list_records[0].seq)
        id_seq_list = [id,seq]
    return id_seq_list

#It parses the primer_barcode info file to a dictionary
def parse_primer_barcode(input_info):
    input_info_dict = {}
    with open(input_info) as f:
        for line in f:
            line = line.rstrip() 
            line_list = line.split("\t")
            input_info_dict[line_list[0]] = [line_list[1],line_list[2],line_list[3],line_list[4]]
    return input_info_dict

#Runs simlord to generate CCS reads
def generate_ccs_reads(input_fasta,output_dir,seq_len,passes_distribution,min_cov,max_cov):
    list_ccs_output_bam = []
    list_intermediate_files_to_remove = []

    output_prefix = os.path.splitext(os.path.basename(input_fasta))[0]
    output_simulated_prefix = os.path.join(output_dir,output_prefix)
    #Get a pseudo randome coverage
    number_reads = random.randint(min_cov, max_cov)
    #simlord --read-reference input_fasta --num-reads 10 --fixed-readlength seq_len -pi 0.0592 -pd 0.0301 -ps 0.0527 myreads
    #Error rate from  Weirather et al. 2017 - F1000Res.
    for i in range(number_reads):
        output_prefix_new = output_prefix + "_" + str(i)
        #cmd = ["simlord","--read-reference",input_fasta,"--num-reads",str(number_reads),"--fixed-readlength",str(seq_len),"-pi","0.0592","-pd","0.0301","-ps","0.0527","--no-sam",output_simulated_prefix]
        number_passes = random.choice(passes_distribution)
        #/project/pbarc_anastrepha_phylo/carlos/bin/pbsim/pbsim3/build/bin/pbsim --strategy templ --method qshmm --qshmm /project/pbarc_anastrepha_phylo/carlos/bin/pbsim/pbsim3/data/QSHMM-RSII.model --template ../fasta_primer_barcode/human_1.1.fasta  --pass-num 10
        #Generate subreads (passes) to create one HiFi read
        cmd = ["pbsim","--strategy","templ","--method","qshmm","--qshmm","/project/pbarc_anastrepha_phylo/carlos/bin/pbsim/pbsim3/data/QSHMM-RSII.model","--template",input_fasta,"--pass-num",str(number_passes),"--prefix", output_prefix_new]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, cwd= output_dir)
        out = p.communicate()
        #Convert SAM to BAM
        output_sam = output_simulated_prefix + "_" + str(i) + ".sam"
        list_intermediate_files_to_remove.append(output_sam)
        output_bam = output_simulated_prefix + "_" + str(i) + ".bam"
        list_intermediate_files_to_remove.append(output_bam)
        cmd = ["samtools","view","-b","-o",output_bam,output_sam]
        p = subprocess.Popen(cmd, cwd= output_dir)
        out = p.communicate()
        #Get the HiFi read
        output_CCS_bam = output_simulated_prefix + "_" + str(i) + "_ccs.bam"
        cmd = ["ccs","-j","1",output_bam,output_CCS_bam]
        p = subprocess.Popen(cmd, cwd= output_dir)
        out = p.communicate()
        list_ccs_output_bam.append(output_CCS_bam)
        list_intermediate_files_to_remove.append(output_CCS_bam)

        #Addional intermediate files
        output_maf = output_simulated_prefix + "_" + str(i) + ".maf"
        list_intermediate_files_to_remove.append(output_maf)
        output_CCS_bam_pbi = output_simulated_prefix + "_" + str(i) + "_ccs.bam.pbi"
        list_intermediate_files_to_remove.append(output_CCS_bam_pbi)
        output_CCS_bam_report = output_simulated_prefix + "_" + str(i) + "_ccs.ccs_report.txt"
        list_intermediate_files_to_remove.append(output_CCS_bam_report)
        output_CCS_bam_json = output_simulated_prefix + "_" + str(i) + "_ccs.zmw_metrics.json.gz"
        list_intermediate_files_to_remove.append(output_CCS_bam_json)

    #Make a file with a list of bam ccs paths
    output_ccs_bam_list = output_simulated_prefix + "_ccs_bam_list"
    with open(output_ccs_bam_list, "w") as ccs_bam_list:
        for bam_ccs_path in list_ccs_output_bam:
            ccs_bam_list.write(bam_ccs_path + "\n")

    #Generates the bam file with the HiFi reads
    output_ccs_bam = output_simulated_prefix + ".bam"
    cmd = ["samtools","merge","-b",output_ccs_bam_list,output_ccs_bam]
    p = subprocess.Popen(cmd, cwd= output_dir)
    out = p.communicate()

    for file in list_intermediate_files_to_remove:
        try:
            os.remove(file)
        except:
            print("Error while deleting a file in post_process function: " + str(file))


    return

#Multitask function of simulate reads
def multitask_get_simulations(list_fasta_barcode_primer,output_dir,seq_len_list,passes_distribution,min_cov,max_cov,number_processors):
    # create pool
    with mp.Pool(min(len(list_fasta_barcode_primer), number_processors)) as p:
        p.starmap(generate_ccs_reads, zip(list_fasta_barcode_primer,itertools.repeat(output_dir),seq_len_list,itertools.repeat(passes_distribution),itertools.repeat(min_cov),itertools.repeat(max_cov)))
    return

#Parse the number of passes file
def parse_passes_distribution(input_passes):
    number_passes = []
    with open(input_passes) as input_passes_file:
        for line in input_passes_file:
            line = line.rstrip()
            number_passes.append(int(line))
    return number_passes

'''
MAIN
'''
def main():
    args = parser.parse_args()
    #Check the flags, get full paths and make directories
    if args.input_dir:
        input_dir = convert_to_full_path(args.input_dir)

    if args.p:
        number_processors = int(args.p)

    if args.min_cov:
        min_cov = int(args.min_cov)

    if args.max_cov:
        max_cov = int(args.max_cov)

    if args.output_dir_fasta:
        output_dir_fasta = convert_to_full_path(args.output_dir_fasta)
        create_dir(output_dir_fasta)
    
    if args.output_dir_simulations:
        output_dir_simulations = convert_to_full_path(args.output_dir_simulations)
        create_dir(output_dir_simulations)

    if args.input_info_primer_barcode:
        input_info_primer_barcode = convert_to_full_path(args.input_info_primer_barcode)
    
    if args.input_passes_dist:
        input_passes_dist = convert_to_full_path(args.input_passes_dist)

    #Parse primer_barcode file
    input_info_dict = parse_primer_barcode(input_info_primer_barcode)
    input_passes_dist_list = parse_passes_distribution(input_passes_dist)

    #Create lists
    list_fasta_barcode_primer = []
    seq_len_list = []
    #Iterate the input dir
    for fasta in glob(os.path.join(input_dir, '*.fasta')):
        id_seq_list = parse_fasta(fasta)
        entry_info = input_info_dict[id_seq_list[0]]
        sequence = id_seq_list[1].upper()
        primer_f = entry_info[0].upper()
        primer_r = str(Seq(entry_info[1].upper()).reverse_complement())
        barcode_f = entry_info[2].upper()
        barcode_r = str(Seq(entry_info[3].upper()).reverse_complement())
        new_seq = barcode_f + primer_f + sequence + primer_r + barcode_r
        output_file = os.path.basename(fasta)
        output_fasta = os.path.join(output_dir_fasta, output_file)
        with open(output_fasta,"w") as out_fasta:
            out_fasta.write(">" + id_seq_list[0] + "\n" + new_seq + "\n")
        list_fasta_barcode_primer.append(output_fasta)
        seq_len_list.append(len(new_seq))

    multitask_get_simulations(list_fasta_barcode_primer,output_dir_simulations,seq_len_list,input_passes_dist_list,min_cov,max_cov,number_processors)
    return

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
    parser = argparse.ArgumentParser(description="This program generates simulated HiFi reads using PBSIM3. The following programs must be available in the path: PBSIM3, samtools, CCS ")
    parser.add_argument('--input_dir', help='Path of the directory containing one or more fasta files to simulate the CCS reads from.')
    parser.add_argument('--input_passes_dist', help='Path of a file with number of passes per HiFi read sampled from a real PacBio run.')
    parser.add_argument('--p', help='Number of processors. CCS may use a considerable amount of memory (10G of memory). For example, it is advisable to use 9 threads for a computer with 100GB RAM.')
    parser.add_argument('--output_dir_fasta', help='Path to the directory to save the fasta files.')
    parser.add_argument('--output_dir_simulations', help='Path to the directory to save the simulated files')
    parser.add_argument('--input_info_primer_barcode', help='Path to file containing the information of the primers and barcode to be added to the sequence. The file must contain the fields: id (exact same as in the fasta file), sequence of primer forward, sequence of primer reverse, sequence of barcode_f, sequence of barcode_r. ')
    parser.add_argument('--min_cov', help='Minimum coverage (HiFi reads) per simulated sequence. To fix the coverage select the same number as in max_cov.')
    parser.add_argument('--max_cov', help='Maximum coverage (HiFi reads) per simulated sequence. To fix the coverage select the same number as in min_cov.')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    #args = parser.parse_args()

    # run main
    main()
