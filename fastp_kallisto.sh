#!/bin/bash

################################################################################

##### USER INPUTS #####
#######################

### Step 1) Define a base directory for pipeline outputs
base_dir=/your/base/directory

### Step 2) Define the two mate pair gzipped FASTQ files and a sample name
mate1=/path/to/file/filename1.fq.gz
mate2=/path/to/file/filename2.fq.gz
sample_name="name_of_sample" # No spaces allowed

### Step 3) Provide the path and filename of a kallisto index of the target
### reference genome
index=/path/to/index/index_name.k31

### Step 4) Set the minimum read length. Reads shorter than "min_len" after
### trimming will be discarded. The default value for fastp is 15 bp.
min_len=15

### Step 5) Set the number of threads to use per sample pair. I have found one
### thread is adequate for the wheat genome. Must be an integer.
threads=1

### Step 6) Ensure your bash implementation (whether on your local machine or 
### a cluster environment) can access a copy of Fastp (>=v0.20.0) and Kallisto
### (>=v0.42.2). You may also convert this script into a batch script for your
### cluster's job handler, e.g. SLURM.


################################################################################

##### SETUP #####
#################

date
echo

### Set up directories for pipeline

out_dir=$base_dir/outputs
mkdir -p $base_dir
mkdir -p $out_dir
mkdir -p $out_dir/fastp_reports

################################################################################

##### RUN SCRIPTS #####
#######################

fastp -i $mate1 \
-I $mate2 \
-o $out_dir/${sample_name}_trimmed_1.fq.gz \
-O $out_dir/${sample_name}_trimmed_2.fq.gz \
-h $out_dir/fastp_reports/${sample_name}_fastp_report.html \
-j $out_dir/fastp_reports/${sample_name}_fastp_report.json \
--detect_adapter_for_pe \
--length_required $min_len \
--thread $threads

echo
echo "Fastp trimming and clean-up complete"
echo

kallisto quant --index $index \
--output-dir $out_dir/$sample_name \
--bootstrap-samples=30 \
--bias \
--threads $threads \
$out_dir/${sample_name}_trimmed_1.fq.gz \
$out_dir/${sample_name}_trimmed_2.fq.gz

echo "Kallisto pseudoalignment complete"
echo
date
