#!/bin/bash
#PBS -N benchmarkz_DD_experiment
#PBS -A UCDV0023
#PBS -l walltime=12:00:00
#PBS -q main
#PBS -l job_priority=regular
#PBS -j oe
#PBS -k eod
#PBS -l select=20:ncpus=128:mpiprocs=128:mem=235GB

export PROSE_REPO_PATH
source ../../scripts/set_env.sh

python3 ${PROSE_REPO_PATH}/scripts/prose_search.py -s setup.ini