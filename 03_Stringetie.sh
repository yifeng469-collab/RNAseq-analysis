#!/bin/bash
#PBS -N stringtie_muscle
#PBS -l nodes=jxlab04:ppn=1,mem=10gb
#PBS -e /home/ZYFLAB/qsublog/
#PBS -o /home/ZYFLAB/qsublog/
#PBS -q cu01
#PBS -t 201-295

# Kill script if any commands fail
set -e

INPUT=/safedisk3/ZYFLAB/projects/F6_RNA_seq/FinalCleandata_sortedbam_muscle/sortedbam
output=/home/ZYFLAB/projects/BackFat_IMF_RNAseq/01_Gene_expression/ALL_Tissues/stringtie
stringTie=/home/zhuyl/software/stringtie-1.3.3b.Linux_x86_64/stringtie
gtf=/home/ZYFLAB/data/RefGenome/susScr11_STAR/Sus_scrofa.Sscrofa11.1.90.gtf
REF=/home/ZYFLAB/data/RefGenome/susScr11_STAR/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa


cd $INPUT

bam=(`ls *.sorted.bam`)

echo "ARRAY JOB: ${bam[$PBS_ARRAYID]}"
NAME=`echo ${bam[$PBS_ARRAYID]} | awk '{split($0,arra,".sorted.bam"); print arra[1]}'`


cd $output

$stringTie ${INPUT}/${bam[$PBS_ARRAYID]} -b ${output} -e -G $gtf -C ${output}/${NAME}_cov_ref.gtf -p 20 -o ${output}/${NAME}
