#!/bin/bash
#PBS -N featureCount
#PBS -l nodes=jxlab03:ppn=1,mem=6gb
#PBS -e /home/ZYFLAB/qsublog/
#PBS -o /home/ZYFLAB/qsublog/
#PBS -q cu01
#PBS -t 184-295

# Kill script if any commands fail
set -e

echo "Job started at:"
date

nprocs=`wc -l < $PBS_NODEFILE`

GTF=/home/ZYFLAB/projects/BackFat_IMF_RNAseq/01_Gene_expression/ALL_Tissues/stringtie/Reformed_merged.gtf
output=/home/ZYFLAB/projects/BackFat_IMF_RNAseq/01_Gene_expression/ALL_Tissues/featureCounts
sortedBAM=/safedisk3/ZYFLAB/projects/F6_RNA_seq/FinalCleandata_sortedbam_muscle/sortedbam
featureCounts=/home/zhuyl/software/subread-1.5.0-p2-Linux-x86_64/bin/featureCounts

cd $sortedBAM
bam=(`ls *.sorted.bam`)

echo "ARRAY JOB: ${bam[$PBS_ARRAYID]}"
NAME=`echo ${bam[$PBS_ARRAYID]} | awk '{split($0,arra,".sorted.bam"); print arra[1]}'`

cd $output
${featureCounts} -T ${nprocs} -F GTF -t exon -g gene_id -s 2 -Q 20 -C -p -D 1000 -a ${GTF} -o ${NAME}.fcounts ${sortedBAM}/${NAME}.sorted.bam

echo "Job ended at:"
date
