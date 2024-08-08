#!/bin/bash
#PBS -N MOM6_timing_job
#PBS -A UCDV0023
#PBS -l walltime=00:20:00
#PBS -q main
#PBS -l job_priority=regular
#PBS -j oe
#PBS -k eod
#PBS -l select=10:ncpus=128:mpiprocs=128

nodes_list=$(cat $PBS_NODEFILE | sort | uniq | cut -d''. -f1)

for node in ${nodes_list}; do
    
    [[ "$(hostname)"] == ${node} ]] && continue

    cmd="cd ${PWD}/.. \\
    && rsync -av --exclude='*prose*' --exclude='__*' --exclude='*.ipynb' benchmark_zstar/ benchmark_zstar_${node} \\
    && cd benchmark_zstar_${node} \\
    && source ../../scripts/set_env.sh \\
    && mpiexec --hosts ${node} -n 128 ../../build/intel/ocean_only/repro/MOM6 >> out_${node}.txt \\
    && mv out_${node}.txt ../benchmark_zstar/ \\
    && cd .. \\
    && rm -r benchmark_zstar_${node}"

    ssh -n ${node} "${cmd}" &
done

echo "** waiting"

wait $(jobs -p)