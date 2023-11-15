# Uauy_Lab_RNAseq_mapping_pipeline
A small pipeline for generating transcriptomic data from raw FASTQ files via Fastp trimming and Kallisto pseudoalignment.
Fastp automatically detects adapters and trims them from the raw reads. Any reads smaller than a specified threshold after trimming are discarded.
Kallisto has been shown to assign reads to the correct wheat homoeolog with high accuracy through the use of nullitetrasomic lines<sup>1</sup>.
###### <sup>1</sup>*Borrill et al., 2016 Plant Physiology*

## Analysis of a single RNA-seq sample
This is a simple bash script that runs Fastp and Kallisto on a pair of files which are the paired-end reads of a single RNA-seq sample.
The script contains a clearly demarked "user inputs" section specifying the inputs required for use with the user's own data. Please alter these as needed before running the script.
Note that this script is for UNIX-based operating systems and for the BASH shell.

## Analysis of many RNA-seq samples
This script has been set up to run on the Norwich Bioscience Institutes (NBI) computing cluster through the job handler SLURM. It sets up an array of jobs, one per sample.
The script contains a clearly demarked "user inputs" section specifying the inputs required for use with the user's own data. Please alter these as needed before running the script.

This script can be modified for use with different insitutions' clusters by altering the following:
1) The job handling "#SBATCH" lines need to be removed and replaced with commands appropriate to your cluster's job handler
2) The "source package XXXXX" lines load in the requisite software (Fastp and Kallisto). These need to be altered for use on other clusters, even if they use SLURM, as the IDs assigned to these packages on other clusters will differ.
Note that this script is for UNIX-based operating systems and for the BASH shell.
