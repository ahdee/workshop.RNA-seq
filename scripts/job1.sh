#!/bin/bash
#SBATCH --job-name=example
#SBATCH --mem=1GB
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --error=job1.error
#SBATCH --output=job1.out
#SBATCH --chdir ./test

module load STAR
module load ngsutilsj/0.4.7
gunzip -c ~/workshop.RNA-seq/seq/test_R1.fastq.gz | head -4000 | grep 0:GTGAAA | wc -l > results.txt
echo "Hello world"
pwd
echo "done!"
sleep 180
