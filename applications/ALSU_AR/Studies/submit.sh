#!/bin/bash

## This script executes 100 instances of the correction chain defined in 'main.m' on the alsacc cluster
## Results are saved in 'results.mat' in the corresponding subfolder created by slurm.

#SBATCH --job-name=ALSU_AR_CHAIN
#SBATCH --array=1-100
#SBATCH --partition=alsacc
#SBATCH --constraint=alsacc_c24
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

#SBATCH --output=out1
#SBATCH --error=err1

module load tools/matlab/r2017b


TASK_ID=`printf %s_%s $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID`;

echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "SLURM_ARRAY_JOB_ID:  $SLURM_ARRAY_JOB_ID"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "SLURM_RESTART_COUNT: $SLURM_RESTART_COUNT";
echo "TASK_ID: $TASK_ID";

#(cd ~/; ~/lmutil stat -a)

MDIR=`pwd`
RUNDIR="./$SLURM_ARRAY_JOB_ID/`printf '%05d' $SLURM_ARRAY_TASK_ID`/"
mkdir -p $RUNDIR
cp -n $MDIR/$1.m $RUNDIR/../
cd $RUNDIR
matlab -nodisplay -nosplash -r "disp(pwd); addpath('../');  $1; save('results.mat','results'); exit;" >out 2>err
