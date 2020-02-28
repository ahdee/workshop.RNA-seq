#!/bin/bash

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
echo "this is what I'm calling it ${outfile}"


gunzip -c ${in1} | head -4000 | grep 0:GTGAAA | wc -l