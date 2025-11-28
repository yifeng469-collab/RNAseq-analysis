#!/bin/bash
#PBS -N counts_geneLength
#PBS -l nodes=jxlab03:ppn=1,mem=1gb
#PBS -e /home/ZYFLAB/qsublog/
#PBS -o /home/ZYFLAB/qsublog/
#PBS -q cu01
#PBS -t 0-554
# Kill script if any commands fail
set -e

echo "job started at:"
date


input=/home/ZYFLAB/projects/BackFat_IMF_RNAseq/01_Gene_expression/ALL_Tissues/featureCounts
output=/home/ZYFLAB/projects/BackFat_IMF_RNAseq/01_Gene_expression/ALL_Tissues/Gene_Counts

cd $input
counts=(`ls *.fcounts`)

echo "ARRAY JOB: ${counts[$PBS_ARRAYID]}"
NAME=`echo ${counts[$PBS_ARRAYID]} | awk '{split($0,arra,"."); print arra[1]}'`

cd $output
awk '{print $1"\t"$7}' ${input}/${NAME}.fcounts | sed '1,2d' > ${NAME}.counts.txt
awk '{print $1"\t"$6}' ${input}/${NAME}.fcounts | sed '1,2d' > ${NAME}.geneLength.txt
