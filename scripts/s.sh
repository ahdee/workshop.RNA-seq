#!/bin/bash
#SBATCH --job-name=example
#SBATCH --mem=70GB
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --error=test2.error
#SBATCH --output=test2.out
#SBATCH --chdir ./test

module load STAR
module load ngsutilsj/0.4.7
gunzip -c ../../seq/test_R1.fastq.gz | head -4000 | grep 0:TTAGGC | wc -l > results.txt
echo "Hello world"
pwd
echo "done!"
sleep 180
