#!/bin/bash
#SBATCH --partition=nbi-medium
#SBATCH --mem 24G
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00
#SBATCH --output=Z/%x_%A_%a.out
#SBATCH --error=Z/%x_%A_%a.err
#SBATCH --array=0-Y:2

################################################################################

##### USER INPUTS #####
#######################

### Step 1) Define a base directory for pipeline outputs
base_dir=/your/base/directory

### Step 2) Define the input directory which contains all your gzipped
### raw FASTQ files, with paired files of the format:
### x_1.fq.gz   x_2.fq.gz
### where 'x' is the sample name.
input_dir=/your/input/directory

### Step 3) Provide the path and filename of a kallisto index of the target
### reference genome
index=/path/to/index/index_name.k31

### Step 4) In the "SBATCH --output" and "SBATCH --error" lines, replace 'Z'
### with the directory you wish to store logs in. This directory must
### already exist.

### Step 5) Set 'Y' in the "SBATCH --array" line such that
### Y = (number of samples * 2) - 1
### e.g. for 21 RNA-seq samples, Y = 41
### This is to account for each sample having a pair of read files

### Step 6) Set the minimum read length. Reads shorter than "min_len" after
### trimming will be discarded. The default value for fastp is 15 bp.
min_len=15

### Step 7) Set the number of threads to use per sample pair. I have found one
### thread is adequate for the wheat genome. Must be an integer.
threads=1

################################################################################

##### SETUP #####
#################

echo date
echo

### Set up generic directories for pipeline
out_dir=$base_dir/outputs
mkdir -p $base_dir
mkdir -p $out_dir
mkdir -p $out_dir/fastp_reports
echo "Directories created"

### Source fastp v0.23.1
source package 7b4623ee-93ca-4c66-b740-361a5e43c0e4
### Source kallisto 0.44.0
source package c1ac247d-c4d1-4747-9817-bf03617f979b
echo "Packages loaded"

### Create shorthand SLURM array variable
i=$SLURM_ARRAY_TASK_ID

### Create array of input file names
readarray -t file_list < <(ls $input_dir)

### Pull the filename for mate 1 for this array task
sample_mate1_filename=${file_list[ $i ]}

### Strip off trailing "_1.fq.gz" to obtain sample name
sample_prefix=${sample_mate1_filename%_1.fq.gz}

echo "Variables set"

################################################################################

##### RUN SCRIPTS #####
#######################

fastp -i $input_dir/${sample_prefix}_1.fq.gz \
-I $input_dir/${sample_prefix}_2.fq.gz \
-o $out_dir/${sample_prefix}_trimmed_1.fq.gz \
-O $out_dir/${sample_prefix}_trimmed_2.fq.gz \
-h $out_dir/fastp_reports/${sample_prefix}_fastp_report.html \
-j $out_dir/fastp_reports/${sample_prefix}_fastp_report.json \
--detect_adapter_for_pe \
--length_required $min_len \
--thread $threads

echo
echo "Fastp trimming and clean-up complete"
echo

kallisto quant --index $index \
--output-dir $out_dir/$sample_prefix \
--bootstrap-samples=30 \
--bias \
--threads $threads \
$out_dir/${sample_prefix}_trimmed_1.fq.gz \
$out_dir/${sample_prefix}_trimmed_2.fq.gz

echo "Kallisto pseudoalignment complete"
echo
echo date
