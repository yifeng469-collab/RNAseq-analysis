#!/bin/bash
#PBS -N sortbam
#PBS -l nodes=1:ppn=2,mem=10gb
#PBS -e /home/ZYFLAB/qsublog/
#PBS -o /home/ZYFLAB/qsublog/
#PBS -q fat
#PBS -t 0-1

# Kill script if any commands fail
set -e

echo "job started at:"
date


bamfile=/home/ZYFLAB/projects/F6_RNAseq/STAR
output=/home/ZYFLAB/projects/F6_RNAseq/sortedbam
samtools=/home/zhangjunjie/bin/samtools-1.3/samtools
TMP=/home/ZYFLAB/projects/F6_RNAseq/TMP
cd $bamfile
bam=(`ls *.bam`)

echo "ARRAY JOB: ${sam[$PBS_ARRAYID]}"
NAME=`echo ${bam[$PBS_ARRAYID]} | awk '{split($0,arra,"."); print arra[1]}'`

echo "sort bam file"
${samtools} sort -T ${TMP}/${NAME}.sorted -o ${output}/${NAME}.sorted.bam ${bamfile}/${NAME}.bam 

echo "create indexed bam file:" 
${samtools} index ${output}/${NAME}.sorted.bam

echo "job finished at:"
date

