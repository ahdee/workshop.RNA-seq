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
echo "this is what I'm calling it ${outfile}"



## run fastqc
fastqc "${in1}" "${in2}" \
--outdir=$output

## run trimm 
java -Xmx15g -jar /asclab/users/leea3/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
"${in1}" \
"${in2}" \
-baseout "${output}${outfile}.fq.gz" \
ILLUMINACLIP:/asclab/users/leea3/bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 \
-threads 6

## run STAR 
STAR --genomeDir /asclab/projects/spcg/resources/indexes/GRCh38.p5/star-index-2.5.3a \
--runThreadN 6 \
--readFilesIn \
"${output}${outfile}_1P.fq.gz" \
"${output}${outfile}_2P.fq.gz" \
--readFilesCommand zcat \
--outFileNamePrefix "${output}${outfile}" \
--outSAMtype BAM \
SortedByCoordinate --quantMode \
GeneCounts --sjdbGTFfile /asclab/data1/backup/projects/spcg/resources/indexes/GRCh38.p5/gencode.v24.annotation.gtf \
--twopassMode Basic

## run aligment QC
ngsutilsj bam-stats --silent --unique \
--gtf /asclab/data1/backup/projects/spcg/resources/indexes/GRCh38.p5/gencode.v24.annotation.gtf \
-o "${output}${outfile}.stat" \
"${output}${outfile}Aligned.sortedByCoord.out.bam"

## create index
samtools index "${output}${outfile}Aligned.sortedByCoord.out.bam"

# if you run this this is what will happen 
# sh workflow1.sh ~/workshop.RNA-seq/seq/ ~/workshop.RNA-seq/seq/testout/ test_R1 tt1