#!/bin/bash

#SBATCH --account=bbkr-delta-cpu
#SBATCH --job-name="run_full_job"
##SBATCH --output="a.out.%j.%N.out"
#SBATCH --partition=cpu
#SBATCH --mem=208G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --cpus-per-task=1  
#SBATCH --no-requeue
#SBATCH -t 24:00:00
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out

export NUMPROC=$SLURM_NTASKS
export EVTSTART=0
export EVTEND=5
export CENT=05
export NSAMPLES=50
export OUTPUT=/work/hdd/bbkr/kpala/convergence/${CENT}/${NSAMPLES}
export TMP=${WORKDIR}/tmp
export LOG=${OUTPUT}.log

# Paths for database and configuration
export DATABASE_PATH=/projects/bbkr/kpala/events.db
export CONFIG_PATH=/projects/bbkr/kpala/configs/pbpb${CENT}_${NSAMPLES}.yml

mkdir -p $TMP
mkdir -p ${WORKDIR}/run_tables

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

# Parallel execution using GNU parallel
seq $EVTSTART $EVTEND | parallel -j${NUMPROC} --workdir ${WORKDIR} --results ${LOG}-parallel \
    python wrapper/main.py {} $DATABASE_PATH $CONFIG_PATH
