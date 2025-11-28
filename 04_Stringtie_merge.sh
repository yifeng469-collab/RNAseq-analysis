#!/bin/bash
#PBS -N stringTie_merge
#PBS -l nodes=jxlab01:ppn=2,mem=8gb
#PBS -e /home/ZYFLAB/qsublog/
#PBS -o /home/ZYFLAB/qsublog/
#PBS -q cu01

dir=/home/ZYFLAB/projects/BackFat_IMF_RNAseq/01_Gene_expression/ALL_Tissues/stringtie
stringTie=/home/zhuyl/software/stringtie-1.3.3b.Linux_x86_64/stringtie
gtf=/home/ZYFLAB/data/RefGenome/susScr11_STAR

cd $dir
find . -type f -name "*.gtf" | xargs ls | awk -F '_' '{print$1}' > mergelist.txt
sed -i "s#./##g" mergelist.txt

${stringTie} --merge -p 5 -G ${gtf}/Sus_scrofa.Sscrofa11.1.90.gtf -o stringtie_merged.gtf mergelist.txt
