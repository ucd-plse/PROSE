#!/bin/bash
#PBS -N MPAS_timing_job
#PBS -A UCDV0023
#PBS -l walltime=00:30:00
#PBS -q main
#PBS -l job_priority=regular
#PBS -j oe
#PBS -k eod
#PBS -l select=10:ncpus=64:mpiprocs=64

nodes_list=$(cat $PBS_NODEFILE | sort | uniq | cut -d''. -f1)

for node in ${nodes_list}; do
    
    [[ "$(hostname)"] == ${node} ]] && continue

    cmd="cd ${PWD}/.. \\
    && rsync -av --exclude='*prose*' --exclude='__*' --exclude='*.ipynb' 240km_uniform/ 240km_uniform_${node} \\
    && cd 240km_uniform_${node} \\
    && source ../../scripts/set_MPAS_env_intel.sh \\
    && make clean \\
    && mpirun --hosts ${node} -n 64 ../../atmosphere_model \\
    && mv log.atmosphere.0000.out ../240km_uniform/out_${node}.txt \\
    && cd .. \\
    && rm -r 240km_uniform_${node}"

    ssh -n ${node} "${cmd}" &
done

echo "** waiting"

wait $(jobs -p)