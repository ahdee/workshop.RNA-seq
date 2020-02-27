#!/bin/bash
#SBATCH --job-name=submit1
#SBATCH --mem=70GB
#SBATCH --ntasks=6
#SBATCH --time=10:00:00
#SBATCH --error=sample1.error
#SBATCH --output=sample2.out
#SBATCH --chdir ./testout

module load STAR/2.5.3a
module load ngsutilsj/0.4.7
module load samtools/1.2

base="$1" # setting up my base directory
output="$2"
inputfile="$3"
outfile="$4"
in1="${base}${inputfile}.fastq.gz"
in2=$(echo $in1 | sed s/"R1"/"R2"/)

# create output dir if it does not exists
mkdir -p $output

# you can AND should define all your variables up here instead of hard coding it. 
#trim.adap="path to stuff"

echo "this is where my seqs are at: $base"
echo "this is where my outputs will be: $output"
echo "this is where my input files are ${in1} ${in2}"
echo "this is we are calling it ${outfile}"