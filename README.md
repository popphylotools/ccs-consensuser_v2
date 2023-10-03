
<p align="left">
<img src="https://github.com/popphylotools/ccs-consensuser_v2/blob/main/ccs-consensuser_logo.jpg" width="700" />
</p>

Created By: Carlos Congrains and Scott Geib. Email: carloscongrains@gmail.com

###################################

CCS-consensuser is a pipeline written in Python 3 that aims to generate consensus of unique sequences (haplotypes) from long-read amplicon sequences obtained by PacBio circular consensus sequencing (CCS). This tool is implemented in two main steps, demultiplexing and consensus.

<p align="center">
<img src="https://github.com/popphylotools/ccs-consensuser_v2/blob/main/ccs-consensuser_v2_workflow.jpg" width="700" />
</p>

REQUIREMENTS AND INSTALLATION:
-----
We developed and tested this script using Python version 3. It requires numpy and Biopython v.1.72 modules. The following programs should be previously installed: 

* mothur v.1.46.1 (https://mothur.org)
* blast v.2.5 (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.5.0/)
* mafft v.7.475 (https://mafft.cbrc.jp/alignment/software/)
* muscle v.3.8.1551 (http://www.drive5.com/muscle/)
* clustalw v.2.1 (http://www.clustal.org/clustal2/)

We strongly recommend to use conda (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to install the dependecies of this program. Once conda is installed in the system, run the following commands to install CCS-consensuser:

```sh
git clone https://github.com/popphylotools/ccs-consensuser_v2.git
cd 
conda env create -f environment.yml
```

USAGE:
-----

Uncompress the test_data.tar.gz to test the program:

```sh
tar -xzvf test_data.tar.gz
```

MODE DEMULTIPLEX:
-----

This mode implements the reorientation of the reads and demultiplexing. A general command line is shown below:

python ccs-consensuser_v2.py --number_processors 10 --mode demultiplex -oligos path_to_oligos_file -i path_to_hifi_reads.fastq -o path_to_output_directory

To test this mode, run the following command line in the place you download this program:

```sh
python ccs-consensuser_v2.py --number_processors 10 --mode demultiplex -i test_data/fastq_oligos/nuclear_loci.fastq -oligos test_data/fastq_oligos/oligos_nuclear -o mode_demultiplex_output 2> mode_demultiplex_output.log
```

The output consists in a set of fastq files with demultiplexed reads. The name files will have the following structure:

basename_input_file_oriented.sample_name.locus_name.fastq

For example:

nuclear_loci_oriented.B_zonata_Nepal_2643.2.PB_o8713_158-707.fastq

MODE CONSENSUS:
-----

This mode aims to obtain a consensus from reads in fastq format from a demultiplexed sample. A general command line is shown below:

python ccs-consensuser_v2.py --mode consensus -l mafft -F sequence_primerF -R sequence_primerR --out_dir path_to_output_directory  -i path_to_hifi_reads.fastq -pm 2 -bm 2 --min_base_score 30 --min_seq_score 20 --sequence_max_mask 10 --min_seq_cov 2 --min_SNP_cov 0 --max_len_delta 5 --alignment_max_amb 0 --consensus_threshold 0.7 --consensus_ignore_mask_char --min_prop_reads 0.2 --mask_char N --indel F 2> path_to_output_log

To test this mode, run the following command line in the place you downloaded this program:

```sh
python ccs-consensuser_v2.py --mode consensus -l mafft -F GCAGTCGAACATGTAGCTGACTCAGGTCACGGGATGGATGTATGTCTGCTG -R TGGATCACTTGTGCAAGCATCACATCGTAGCAAATACCAAAACCCCATAGCCAT --out_dir mode_consensus_output -i  test_data/mode_consensus/B_tau_Laos_330.PB_o5188_853-1541.fastq -pm 2 -bm 2 --min_base_score 30 --min_seq_score 20 --sequence_max_mask 10 --min_seq_cov 2 --min_SNP_cov 0 --max_len_delta 5 --alignment_max_amb 0 --consensus_threshold 0.7 --consensus_ignore_mask_char --min_prop_reads 0.2 --mask_char N --indel F 2> mode_consensus_output.log
```

The output consists in a set of fasta files and a zip file (--clean_up zip). For instance, the script will generate the following files after running the command line indicated above:

B_tau_Laos_330.PB_o5188_853-1541_1.3.consensus.fasta
B_tau_Laos_330.PB_o5188_853-1541_18.2.consensus.fasta
B_tau_Laos_330.PB_o5188_853-1541_7.6.consensus.fasta
mode_consensus_output.1.minor_freq_hap.consensus.fasta
mode_consensus_output.2.consensus.fasta
mode_consensus_output.2.consensus.aln.fasta
mode_consensus_output_intermediate_files.zip

The three first files are the consensus of each haplotype found in the sample. The filenames inform the sample name (B_tau_Laos_330), marker name (PB_o5188_853-1541) and coverage (3,2 and 6).
The fourth file contains the haplotype(s) that did not pass the filter --min_prop_reads, which in this case was 0.2. The file name inficates the number of haplotypes in this file, which was only 1 (B_tau_Laos_330.PB_o5188_853-1541_18.2.consensus.fasta).
The mode_consensus_output.2.consensus.fasta is the final result containing the sequences of the haplotypes found in this sample.
In case of more than one haplotype, the program also produces an alignment, stored in mode_consensus_output.2.consensus.aln.fasta.
The zip file include all the intermediate files created for this analysis. If the user does not require these files use --clean_up remove or if the user would like to check those files use --clean_up keep.

MODE CONSENSUS_BATCH:
-----
The consensus batch mode is similar to consensus mode, but it performs this task in a parallel way. The user should prepare a file containing the primers, path of the PacBio Hifi fastq files and path of the output directory. The information of each sample is added in separate lines in a file named as input_list. An example of the input_list file is displayed below:

```sh
primer  GCAGTCGAACATGTAGCTGACTCAGGTCACACACGCCATTTTACTATTAAAACGG TGGATCACTTGTGCAAGCATCACATCGTAGGCTTTCATTCCCTGCCATGG     test_data/mode_consensus_batch/PB_o5503_129-904/B_tau_Laos_330/nuclear_loci_example_oriented.B_tau_Laos_330.PB_o5503_129-904.fastq     test_data/mode_consensus_batch/PB_o5503_129-904/B_tau_Laos_330
primer  GCAGTCGAACATGTAGCTGACTCAGGTCACACACGCCATTTTACTATTAAAACGG TGGATCACTTGTGCAAGCATCACATCGTAGGCTTTCATTCCCTGCCATGG     test_data/mode_consensus_batch/PB_o5503_129-904/B_tau_Thailand_354.1/nuclear_loci_example_oriented.B_tau_Thailand_354.1.PB_o5503_129-904.fastq test_data/mode_consensus_batch/PB_o5503_129-904/B_tau_Thailand_354.1
```

A general command line is shown below:

python ccs-consensuser_v2.py --number_processors 10 --mode consensus_batch --clean_up zip --in_file_list path_to_input_list -pm 2 -bm 2 --min_base_score 30 --min_seq_score 20 --sequence_max_mask 10 --min_seq_cov 2 --min_SNP_cov 0 --max_len_delta 5 --alignment_max_amb 0 --consensus_threshold 0.7 --consensus_ignore_mask_char --min_prop_reads 0.2 --mask_char N --indel F 2> path_to_output_log

To test this mode, run the following command line in the place you downloaded this program:

```sh
python ccs-consensuser_v2.py --number_processors 10 --mode consensus_batch --clean_up zip --in_file_list test_data/mode_consensus_batch/consensus_batch_input -pm 2 -bm 2 --min_base_score 30 --min_seq_score 20 --sequence_max_mask 10 --min_seq_cov 2 --min_SNP_cov 0 --max_len_delta 5 --alignment_max_amb 0 --consensus_threshold 0.7 --consensus_ignore_mask_char --min_prop_reads 0.2 --mask_char N --indel F 2> consensus_batch_output.log
```

The output of each sample will be stored in the directory provided in the fifth field of the input_list file. For more information about the output see MODE CONSENSUS output above.

MODE ALL:
-----
This mode runs the complete pipeline, which is reorientation of the reads, demultiplexing and obtain the consensus of each sample. A general command line is shown below:

python ccs-consensuser_v2.py --number_processors 10 --mode all -l mafft -i path_to_hifi_reads.fastq -oligos path_to_oligos_file -o path_to_output_directory -pm 2 -bm 2 --min_base_score 30 --min_seq_score 20 --sequence_max_mask 10 --min_seq_cov 2 --min_SNP_cov 0 --max_len_delta 5 --alignment_max_amb 0 --consensus_threshold 0.7 --consensus_ignore_mask_char --min_prop_reads 0.2 --mask_char N --indel F 2> path_to_output_log

To test this mode, run the following command line in the place you downloaded this program:

```sh
python ccs-consensuser_v2.py --number_processors 10 --mode all -l mafft -i test_data/fastq_oligos/nuclear_loci.fastq -oligos test_data/fastq_oligos/oligos_nuclear -o mode_all_output -pm 2 -bm 2 --min_base_score 30 --min_seq_score 20 --sequence_max_mask 10 --min_seq_cov 2 --min_SNP_cov 0 --max_len_delta 5 --alignment_max_amb 0 --consensus_threshold 0.7 --consensus_ignore_mask_char --min_prop_reads 0.2 --mask_char N --indel F 2> mode_all_output.log
```

This mode will make a directory for each locus contining one directory per sample. For more information about the output see MODE CONSENSUS output above.

OUTPUT:
-----

The output files are organized in directories, one per sample and one referenced to as final_consensus_sequences. When multiple markers are used (primer pairs), the program generates one directory per marker. Example of the output tree structure:

```sh
Output_directory/
├─ marker1/
│  ├─ final_consensus_sequences
│  ├─ sample01
│  ├─ sample02
├─ marker2/
│  ├─ final_consensus_sequences
│  ├─ sample01
│  ├─ sample02
├─ primers.fasta
├─ input_prefix_oriented.fastq (Input of mothur. Contain the oriented reads in fastq format, used as input for mother.)
├─ input_prefix_oriented.group (Output of mothur. Contain the information of the demuliplex step)
├─ input_prefix_oriented.scrap.fastq (Output of mothur. Contain the information of the demuliplex step)
├─ mothur.1.logfile (Output of mothur. Log file)
```

Content of the final_consensus_sequences directory:
- consensus_sequences.fasta - This file includes all the consensus sequences.

Content of the directory of one sample. The hypothetical example describes the directory of the sample human_180:
- human_180.2.consensus.fasta - Prefix has the form SampleID, "." and number of consensus sequences.
- human_180.2.consensus.aln.fasta - When more than one sequence the program aligns the sequences of the human_180.2.consensus.fasta.
- human_HLA_CCS_rq_099_new_header_oriented.human_180_1.13.consensus.fasta - Independent file for a consensus sequence. Prefix has the form: name of the input file, ".", sampleID (human_180), "_", haplotypeID (1), "." and the coverage (13).
- human_HLA_CCS_rq_099_new_header_oriented.human_180_2.32.consensus.fasta - Independent file for a consensus sequence. Prefix has the form: name of the input file, ".", sampleID (human_180), "_", haplotypeID (2), "." and the coverage (32).
- human_180_intermediate_files.zip - A zip compressed file with all the intermediate files.

OPTIONS:
-----



Citation:
-----
Congrains et al 2022.

This script is in the public domain in the United States per 17 U.S.C. § 105
