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

base="${HOME}/workshop.RNA-seq/seq/" # setting up my base directory
output="${HOME}/workshop.RNA-seq/seq/testout/"

# create output dir if it does not exists
mkdir -p $output

# you can AND should define all your variables up here instead of hard coding it. 
#trim.adap="path to stuff"

echo "this is where my seqs are at: $base"
echo "this is where my outputs will be: $output"


## run fastqc
fastqc "${base}test_R1.fastq.gz" "${base}test_R2.fastq.gz" \
--outdir=$output

## run trimm 
java -Xmx15g -jar /asclab/users/leea3/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
"${base}test_R1.fastq.gz" \
"${base}test_R2.fastq.gz" \
-baseout "${output}test_s1.fq.gz" \
ILLUMINACLIP:/asclab/users/leea3/bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 \
-threads 6

## run STAR 
STAR --genomeDir /asclab/projects/spcg/resources/indexes/GRCh38.p5/star-index-2.5.3a \
--runThreadN 6 \
--readFilesIn \
"${output}test_s1_1P.fq.gz" \
"${output}test_s1_2P.fq.gz" \
--readFilesCommand zcat \
--outFileNamePrefix "${output}testS1" \
--outSAMtype BAM \
SortedByCoordinate --quantMode \
GeneCounts --sjdbGTFfile /asclab/data1/backup/projects/spcg/resources/indexes/GRCh38.p5/gencode.v24.annotation.gtf \
--twopassMode Basic

## run aligment QC
ngsutilsj bam-stats --silent --unique \
--gtf /asclab/data1/backup/projects/spcg/resources/indexes/GRCh38.p5/gencode.v24.annotation.gtf \
-o "${output}testS1.stat" \
"${output}testS1Aligned.sortedByCoord.out.bam"

## create index
samtools index "${output}testS1Aligned.sortedByCoord.out.bam"