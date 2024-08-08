#!/bin/bash
#PBS -N funarc_brute_force
#PBS -A UCDV0023
#PBS -l walltime=2:00:00
#PBS -q regular
#PBS -j oe
#PBS -k eod
#PBS -l select=1:ncpus=36:mem=45GB
#PBS -l inception=login

make reset
python3 ${PROSE_REPO_PATH}/scripts/prose_brute_force_search.py -s setup.ini