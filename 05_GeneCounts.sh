#!/bin/bash
#SBATCH -J GeneCounts
#SBATCH -N 1 -n 3
#SBATCH -e /home/ZYFLAB/qsublog/job-%j_%a.err
#SBATCH -o /home/ZYFLAB/qsublog/job-%j_%a.log
#SBATCH -p jxlab
#SBATCH --mem=10Gb
#SBATCH -a 0-1


#source /opt/software/anaconda2/bin/activate
echo "Job Start at `date`"

#####enviroments setting
PBS_ARRAYID=${SLURM_ARRAY_TASK_ID}
array_jobid=${SLURM_ARRAY_JOB_ID}
jobid=${SLURM_JOB_ID}
job_name=${SLURM_JOB_NAME}
workdir=${SLURM_SUBMIT_DIR}
pid=${SLURM_TASK_PID}
nprocs=${SLURM_JOB_CPUS_PER_NODE}
tmpdir=/tmpdisk
##########software path setting#############
#STAR=/safedisk2/ZYFLAB/cp_ZYFLAB_2021_3_9/bin/STAR/STAR-2.5.3a/bin/Linux_x86_64_static/STAR
#DIR=/home/ZYFLAB/Mouse_RNAseq/RefGenome
input=/home/ZYFLAB/Mouse_RNAseq/FeatureCounts
output=/home/ZYFLAB/Mouse_RNAseq/Gene_Counts
#############setting input and output directory##########
#DIR=/home/ZYFLAB/Mouse_RNAseq/RefGenome
#OUTPUT=/home/ZYFLAB/Mouse_RNAseq/RefGenome



#############working directory setting###########
uid="ZYFLAB" #user id
ls_date=`date +m%d%H%M%S` #Set Random date and time
Work_dir=${uid}_${ls_date}_${PBS_ARRAYID}  #
cd ${tmpdir}
echo "the temprory directory is $tmpdir"
mkdir ${Work_dir}
cd ${tmpdir}/${Work_dir}


#####copy data to node computers ###################
#fa=(`ls ${INPUT}/*.fa`)
#echo "current file is ${fa[$PBS_ARRAYID]}"

#cp ${INPUT}/${Samples[$PBS_ARRAYID]}/*.gtz ./.

##################################################################
########## YOUR RUNNING SCRIPTS                      #############
##################################################################
cd $input
counts=(`ls *.fcounts`)

echo "ARRAY JOB: ${counts[$PBS_ARRAYID]}"
NAME=`echo ${counts[$PBS_ARRAYID]} | awk '{split($0,arra,"."); print arra[1]}'`

cd $output
awk '{print $1"\t"$7}' ${input}/${NAME}.fcounts | sed '1,2d' > ${NAME}.counts.txt
awk '{print $1"\t"$6}' ${input}/${NAME}.fcounts | sed '1,2d' > ${NAME}.geneLength.txt

##--outFilterMismatchNoverReadLmax 0.07 allows 10 mismatches in a 150bp read

#########CP BACK AND SAVE YOUR RESULTS #####################
#cp *  ${OUTPUT}/.

######REMOVE tmp directory #################################

#rm -rf ${tmpdir}/${Work_dir}
#conda deactivate
echo "Job finished at:" `date`

