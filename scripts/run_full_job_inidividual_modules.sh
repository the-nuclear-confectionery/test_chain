#!/bin/bash

#SBATCH --account=bbkr-delta-cpu
#SBATCH --job-name="run_full_job"
##SBATCH --output="a.out.%j.%N.out"
#SBATCH --partition=cpu
#SBATCH --mem=208G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=1  
#SBATCH --no-requeue
#SBATCH -t 08:00:00
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out

export NUMPROC=$SLURM_NTASKS
export EVTSTART=0
export EVTEND=9
export OUTPUT=${WORKDIR}/results
export TMP=${WORKDIR}/tmp
export LOG=${OUTPUT}.log

mkdir -p $OUTPUT
mkdir -p $TMP

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

mkdir -p ${WORKDIR}/run_tables

# Adjusted parallel job submission
seq $EVTSTART $EVTEND | parallel -j${NUMPROC} --workdir ${WORKDIR} --results ${LOG}-parallel \
    ./scripts/run_ccake_chain.sh {}
