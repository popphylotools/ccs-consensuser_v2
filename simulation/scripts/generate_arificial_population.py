from Bio import SeqIO
import os,sys,argparse
import csv
import random
import copy

"""
Introduce variants (SNPs or indels) into a reference sequence in fasta format. This script generates random combination of variants and produces two fasta files per individual (one per chromosome)

Table of variants look like (position, reference, alternate)

3    A   T,C
10    C   G
12  -   TCC
15  TGG -
"""

'''
FUNCTIONS
'''

#Pasrse the table of variants
def parse_csv_table_variants(table_variants):
    # opening the CSV file
    with open(table_variants, mode ='r') as file:
    # reading the CSV file
        table_variants_list = list(csv.reader(file, delimiter="\t"))
    return table_variants_list

#Select random variants based on a maximum value
def select_random_variants(table_variants_list,max_number_variants):
    #Number of variants
    number_variants = random.randrange(1,max_number_variants+1)
    #Select random variants
    list_pos_sample_variants = []
    sample_variants = []
    while len(sample_variants) < number_variants:
        sample_variant = random.sample(table_variants_list,1)
        if sample_variant[0] not in sample_variants and sample_variant[0][0] not in list_pos_sample_variants:
            sample_variants.append(sample_variant[0])
            list_pos_sample_variants.append(sample_variant[0][0])

    #sample_variants = random.sample(table_variants_list,number_variants)
    return sample_variants

#Convert the fasta file (reference) into a list (id and sequence)
def parse_fasta(input_fasta):
    list_records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        list_records.append(record)

    if len(list_records) == 1:
        id = list_records[0].id
        seq = str(list_records[0].seq)
        id_seq_list = [id,seq]
    return id_seq_list

def process_fasta(id_seq_list,table_variants_list,max_number_variants,number_ind_pop,output_DIR,prefix):
    #Convert sequence string into a list
    initial_allele_1 = [i for i in id_seq_list[1].upper()]
    initial_allele_2 = [i for i in id_seq_list[1].upper()]

    #Generates the population
    for ind in range(number_ind_pop):
        new_initial_allele_1 = copy.deepcopy(initial_allele_1)
        new_initial_allele_2 = copy.deepcopy(initial_allele_2)
        diploid = [new_initial_allele_1,new_initial_allele_2]
        sample_variants = select_random_variants(table_variants_list,max_number_variants)
        for variant in sample_variants:
            pos = int(variant[0])
            #If multiple alternate alleles, randomly select one
            alt = random.sample(list(variant[2].split(",")),1)[0]
            #Decide if the individual will be homozygous or heterozygous
            state = random.sample(["hom","het"],1)
            if state == ["hom"]:
                diploid[0][pos] = alt
                diploid[1][pos] = alt
            else:
                diploid[random.randrange(2)][pos] = alt

        filename_alelle_1 = os.path.join(output_DIR, prefix + "_" + str(ind) + ".1.fasta")
        filename_alelle_2 = os.path.join(output_DIR, prefix + "_" + str(ind) + ".2.fasta")
        #save the files in fasta format
        with open(filename_alelle_1,"w") as file:
            sequence_str_haplotype1 = "".join(diploid[0])
            new_sequence_str_haplotype1 = ""
            for nucleotide in sequence_str_haplotype1:
                if nucleotide == "-":
                    nucleotide = ""
                new_sequence_str_haplotype1 = new_sequence_str_haplotype1 + nucleotide
            file.write(">" + prefix + "_" + str(ind) + ".1" +"\n" + new_sequence_str_haplotype1 + "\n")
        with open(filename_alelle_2,"w") as file:
            sequence_str_haplotype2 = "".join(diploid[1])
            new_sequence_str_haplotype2 = ""
            for nucleotide in sequence_str_haplotype2:
                if nucleotide == "-":
                    nucleotide = ""
                new_sequence_str_haplotype2 = new_sequence_str_haplotype2 + nucleotide
            file.write(">" + prefix + "_" + str(ind) + ".2" +"\n" + new_sequence_str_haplotype2 + "\n")
    return

#Convert a relative path to full path
def convert_to_full_path(path):
    absolute_path = os.path.abspath(path)
    return absolute_path

#Create a directory
def create_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    return dir_path

#Main function
def main(_args):
    #Convert to full paths
    if _args.input_fasta:
        input_fasta = convert_to_full_path(_args.input_fasta)
    if _args.table_variants:
        table_variants = convert_to_full_path(_args.table_variants)
    if _args.output_DIR:
        output_DIR = convert_to_full_path(_args.output_DIR)
    if _args.prefix:
        prefix = _args.prefix
    if _args.max_number_variants:
        max_number_variants = _args.max_number_variants
    if _args.number_ind_pop:
        number_ind_pop = _args.number_ind_pop

    create_dir(output_DIR)
    table_variants_list = parse_csv_table_variants(table_variants)
    id_seq_list = parse_fasta(input_fasta)
    process_fasta(id_seq_list,table_variants_list,max_number_variants,number_ind_pop,output_DIR,prefix)
    return

"""
MAIN
"""

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
    parser.add_argument('--input_fasta', help='Path to the reference sequence in fasta format.')
    parser.add_argument('--sample_prefix', help='Prefix to be used to save the fasta files.')
    parser.add_argument('--table_variants', help='Path to the table of variants containing the following tab delimited fields: position, reference allele, alternate allele(s). Use the comma to separate multiple alternate alleles. This variants can be either SNPs or indels.')
    parser.add_argument('--output_DIR', help='Path to the output directory.')
    parser.add_argument('--max_number_variants', type=int, help='Maximum number of variants to be added into the reference to generate the artificial haplotype. Must be an integer.')
    parser.add_argument('--number_ind_pop', type=int, help='Number of diploid individuals to be generated as part of the artificial population. Must be an integer.')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # run main
    try:
        exit(main(args))
    except Exception as e:
        print("Exception in main(): {}".format(e))
        exit(1)

