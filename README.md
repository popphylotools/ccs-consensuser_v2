
<p align="left">
<img src="https://github.com/popphylotools/ccs-consensuser_v2/blob/main/ccs-consensuser_v2_logo.jpg" width="700" />
</p>


######################################

CCS-consensuser is a pipeline written in Python 3 that aims to generate consensus of unique sequences (haplotypes) from long-read amplicon sequences obtained by PacBio circular consensus sequencing (CCS). This tool is implemented in two main steps, demultiplexing and consensus.

<p align="center">
<img src="https://github.com/popphylotools/ccs-consensuser_v2/blob/main/ccs-consensuser_v2_workflow.jpg" width="700" />
</p>

REQUIREMENTS AND INSTALLATION:
-----
We developed and tested this script using Python version 3. The following dependencies are also required:

Python modules:

* numpy
* biopython v.1.83

Programs with the corresponding versions:

* mothur v.1.46.1 (https://mothur.org)
* blast v.2.16 (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/)
* mafft v.7.490 (https://mafft.cbrc.jp/alignment/software/)
* muscle v.3.8.31 (http://www.drive5.com/muscle/)
* clustalw v.2.1 (http://www.clustal.org/clustal2/)

Installation methods:
---
1. This method requires conda (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) installed in system. The program was tested using conda 24.11.0. Once conda is properly installed, run the following commands to install dependencies of CCS-consensuser:

```sh
#Download the program folder
git clone https://github.com/popphylotools/ccs-consensuser_v2.git
cd ccs-consensuser_v2
#Linux OS (mothur v.1.46.1)
conda env create -f environment_LinuxOS.yml
#macOS. mothur v.1.46.1 is not available in conda for macOS. mothur 1.44.11 will be installed instead.
conda env create -f environment_macOS.yml
```
2. This method requires docker (https://docs.docker.com/engine/install/) installed in the system.

```sh
#Download the container
docker pull carloscongrains/ccs-consensuser
#Run the container
docker run --rm -v`pwd`:`pwd` carloscongrains/ccs-consensuser ccs-consensuser_v2.1.py
```

3. This method requires singularity or apptainer (https://apptainer.org/admin-docs/master/installation.html) installed in the system.

```sh
#Download the container and create the sif file.
singularity pull docker://carloscongrains/ccs-consensuser:2.1
#Run the container. Use the path of the sif file created on the previous step to run the program:
singularity exec ccs-consensuser_2.1.sif ccs-consensuser_v2.1.py
```

4. The users may opt to install the dependencies manually. It is important to install the specified versions of the dependencies. This method may be time consuming.

USAGE:
-----

Uncompress the test_data.tar.gz to test the program:

```sh
tar -xzvf test_data.tar.gz
```

MODE DEMULTIPLEX:
-----

This mode implements the reorientation of the reads and demultiplexing. A general command line is shown below:

python ccs-consensuser_v2.1.py --number_processors 10 --mode demultiplex -oligos path_to_oligos_file -i path_to_hifi_reads.fastq -o path_to_output_directory

To test this mode, run the following command line in the location you downloaded the program and extracted the test_data.tar.gz:

```sh
python ccs-consensuser_v2.1.py --number_processors 10 --mode demultiplex -i test_data/fastq_oligos/nuclear_loci.fastq -oligos test_data/fastq_oligos/oligos_nuclear -o mode_demultiplex_output 2> mode_demultiplex_output.log
```

The output consists in a set of fastq files with demultiplexed reads. The name files will have the following structure:

basename_input_file_oriented.sample_name.locus_name.fastq

For example:

nuclear_loci_oriented.B_zonata_Nepal_2643.2.PB_o8713_158-707.fastq

Additionally, the program creates a summary file refered to "summary_results". This is a tab delimited file with the following fields: MarkerID, SampleID, and Demultiplex_coverage. The summary_results file can be imported in a spreadsheet program such as Microsoft Excel.

MODE CONSENSUS:
-----

This mode aims to obtain a consensus from reads in fastq format from a demultiplexed sample. A general command line is shown below:

python ccs-consensuser_v2.1.py --mode consensus -l mafft -F sequence_primerF -R sequence_primerR --out_dir path_to_output_directory  -i path_to_hifi_reads.fastq -pm 2 -bm 2 --min_base_score 30 --min_seq_score 20 --sequence_max_mask 10 --min_seq_cov 2 --min_SNP_cov 0 --max_len_delta 5 --alignment_max_amb 0 --consensus_ignore_mask_char --min_prop_reads 0.2 --mask_char N --indel M 2> path_to_output_log

To test this mode, run the following command line in the place you downloaded this program:

```sh
python ccs-consensuser_v2.1.py --mode consensus -l mafft -F GCAGTCGAACATGTAGCTGACTCAGGTCACGGGATGGATGTATGTCTGCTG -R TGGATCACTTGTGCAAGCATCACATCGTAGCAAATACCAAAACCCCATAGCCAT --out_dir mode_consensus_output -i  test_data/mode_consensus/B_tau_Laos_330.PB_o5188_853-1541.fastq -pm 2 -bm 2 --min_base_score 30 --min_seq_score 20 --sequence_max_mask 10 --min_seq_cov 2 --min_SNP_cov 0 --max_len_delta 5 --alignment_max_amb 0 --consensus_ignore_mask_char --min_prop_reads 0.2 --mask_char N --indel M 2> mode_consensus_output.log
```

The output consists of a set of fasta files and a zip file (--clean_up zip). For instance, the script will generate the following files after running the command line indicated above:

```sh
B_tau_Laos_330.PB_o5188_853-1541_1.4.consensus.fasta
B_tau_Laos_330.PB_o5188_853-1541_17.2.consensus.fasta
B_tau_Laos_330.PB_o5188_853-1541_7.6.consensus.fasta
mode_consensus_output.1.minor_freq_hap.consensus.fasta
mode_consensus_output.2.consensus.fasta
mode_consensus_output.2.consensus.aln.fasta
mode_consensus_output_intermediate_files.zip
```

The three first files are the consensus of each haplotype found in the sample. The filenames inform the sample name (B_tau_Laos_330), marker name (PB_o5188_853-1541) and coverage, which for this example are 3, 2 and 6, respectively.
The fourth file contains the haplotype(s) that did not pass the filter --min_prop_reads, which in this case was set to 0.2. The file name indicates the number of haplotypes included, which was only 1 (B_tau_Laos_330.PB_o5188_853-1541_18.2.consensus.fasta).
The mode_consensus_output.2.consensus.fasta is the final result for this sample and contains the sequences of the haplotypes in fasta format.
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

python ccs-consensuser_v2.1.py --number_processors 10 --mode consensus_batch --clean_up zip --in_file_list path_to_input_list -pm 2 -bm 2 --min_base_score 30 --min_seq_score 20 --sequence_max_mask 10 --min_seq_cov 2 --min_SNP_cov 0 --max_len_delta 5 --alignment_max_amb 0 --consensus_ignore_mask_char --min_prop_reads 0.2 --mask_char N --indel M 2> path_to_output_log

To test this mode, run the following command line in the place you downloaded this program:

```sh
python ccs-consensuser_v2.1.py --number_processors 10 --mode consensus_batch --clean_up zip --in_file_list test_data/mode_consensus_batch/consensus_batch_input -pm 2 -bm 2 --min_base_score 30 --min_seq_score 20 --sequence_max_mask 10 --min_seq_cov 2 --min_SNP_cov 0 --max_len_delta 5 --alignment_max_amb 0 --consensus_ignore_mask_char --min_prop_reads 0.2 --mask_char N --indel M 2> consensus_batch_output.log
```

The output of each sample will be stored in the directory provided in the fifth field of the input_list file. For more information about the output see MODE CONSENSUS output above.

MODE ALL:
-----
This mode runs the complete pipeline, which is reorientation of the reads, demultiplexing and obtain the consensus of each sample. A general command line is shown below:

python ccs-consensuser_v2.1.py --number_processors 10 --mode all -l mafft -i path_to_hifi_reads.fastq -oligos path_to_oligos_file -o path_to_output_directory -pm 2 -bm 2 --min_base_score 30 --min_seq_score 20 --sequence_max_mask 10 --min_seq_cov 2 --min_SNP_cov 0 --max_len_delta 5 --alignment_max_amb 0 --consensus_ignore_mask_char --min_prop_reads 0.2 --mask_char N --indel M 2> path_to_output_log

To test this mode, run the following command line in the place you downloaded this program:

```sh
python ccs-consensuser_v2.1.py --number_processors 10 --mode all -l mafft -i test_data/fastq_oligos/nuclear_loci.fastq -oligos test_data/fastq_oligos/oligos_nuclear -o mode_all_output -pm 2 -bm 2 --min_base_score 30 --min_seq_score 20 --sequence_max_mask 10 --min_seq_cov 2 --min_SNP_cov 0 --max_len_delta 5 --alignment_max_amb 0 --consensus_ignore_mask_char --min_prop_reads 0.2 --mask_char N --indel M 2> mode_all_output.log
```

This mode will make a directory for each locus containing one directory per sample. For more information about the output see MODE CONSENSUS output above.

To facilitates subsequent analysis, the program generates a fasta file with all consensus sequences per marker. The fasta file of the respective marker is located in subfolders of final_consensus_sequences. The sequences are named based on the information provided in the oligos file and an index (haplotype code). After running the testing command, the program generates two files including the consensus sequences: final_consensus_sequences/PB_o5503_129-904/PB_o5503_129-904.fasta and final_consensus_sequences/PB_o8713_158-707/PB_o8713_158-707.fasta.

Addionally, the program generates a summary file refered to as “summary_results”. This file includes the following fields: Marker identifier, sample identifier, number of demultiplexed reads, number of haplotypes, haplotype identifier, haplotype coverage, number of missing nucleotides (Ns), and haplotype length. The summary_results file can be imported in a spreadsheet program such as Microsoft Excel.

OPTIONS:
-----

<p align="center">
<img src="" width="800" />
</p>

Credits:
-----
[Carlos Congrains](https://github.com/carloscongrains): Contributor  
[Forest Bremer](https://github.com/fbremer): Legacy contributor  
[Julian Dupuis](https://github.com/jrdupuis): Contributor  
[Scott Geib](https://github.com/scottgeib): Contributor  

Citation:
-----
Congrains et al 2024.

This script is in the public domain in the United States per 17 U.S.C. § 105
