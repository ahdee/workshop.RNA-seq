#!/bin/bash
fastqc="fastqc FILE1 FILE2"

for i in 1 2 3 4 5
do
   temp=$(echo $fastqc | sed s/FILE1/NEW_${i}_R1.tar.gz/)
   temp=$(echo $temp | sed s/FILE2/NEW_${i}_R2.tar.gz/)
   echo $temp
done
