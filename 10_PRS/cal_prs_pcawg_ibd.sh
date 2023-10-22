#!/bin/bash
#SBATCH -o cal_pcawg_ibd_prs.out
#SBATCH -J cal_pcawg_ibd_prs
#SBATCH --mem=32G
#SBATCH -c 1
#SBATCH -p short
#SBATCH -t 12:00:00
#SBATCH -a 1-22

cd /n/data1/hms/dbmi/zaklab/csy/NEER_germline_project/10_PRS/DATA/PRS

ln -s /n/app/openblas/0.2.19/lib/liblapack.so liblapack.so.3
export LD_LIBRARY_PATH=$PWD
module load gcc/6.2.0
module load bcftools
module load htslib
module load R/3.6.1


echo "chr$SLURM_ARRAY_TASK_ID"
env SGE_TASK_ID=$SLURM_ARRAY_TASK_ID $*

/home/sc377/nps/sge/nps_score.bgen.job.new ibd.npsdat/ IBD_sung_pcawg/pcawg_bgen pcawg

