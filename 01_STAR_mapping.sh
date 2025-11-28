#!/bin/bash
#PBS -N STAR
#PBS -l nodes=1:ppn=3,mem=10gb
#PBS -e /home/ZYFLAB/qsublog/
#PBS -o /home/ZYFLAB/qsublog/
#PBS -q fat
#PBS -t 0-1


# Kill script if any commands fail
set -e

echo "job started at:"
date

export PATH=/home/zhangjunjie/bin/samtools-1.3:$PATH

REF=/home/ZYFLAB/data/RefGenome/susScr11_STAR
INPUT=/home/ZYFLAB/projects/F6_RNAseq/CleanData
OUTPUT=/home/ZYFLAB/projects/F6_RNAseq/STAR
STAR=/home/ZYFLAB/bin/STAR/STAR-2.5.3a/bin/Linux_x86_64_static/STAR
nprocs=`wc -l < $PBS_NODEFILE`

cd $INPUT
Fastq=(`ls *_R1.fq.gz`)

echo "ARRAY JOB: ${Fastq[$PBS_ARRAYID]}"
NAME=`echo ${Fastq[$PBS_ARRAYID]} | awk '{split($0,arra,"_R1.fq.gz"); print arra[1]}'`

${STAR} --runMode alignReads \
        --runThreadN $nprocs \
        --genomeDir ${REF} \
        --readFilesIn ${INPUT}/${NAME}_R1.fq.gz ${INPUT}/${NAME}_R2.fq.gz \
        --readFilesCommand zcat \
        --clip3pNbases 0 20 \
        --outSJfilterOverhangMin 30 16 16 16 \
        --outSJfilterCountUniqueMin 4 2 2 2 \
        --alignSJoverhangMin 6 \
        --outFilterType BySJout \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --alignIntronMax 500000 \
        --alignMatesGapMax 500000 \
        --outFilterMismatchNoverReadLmax 0.07 \
        --outFileNamePrefix ${OUTPUT}/${NAME} \
        --outStd SAM --outSAMmode Full | samtools view -bS - > ${OUTPUT}/${NAME}.bam

####for trimmed fastq, don't use "--clip3pNbases 0 20" parameter again

echo "job finished at:"
date 

##--outFilterMismatchNoverReadLmax 0.07 allows 10 mismatches in a 150bp read


