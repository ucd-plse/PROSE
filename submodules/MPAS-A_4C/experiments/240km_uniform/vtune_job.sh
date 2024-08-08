#!/bin/bash
#PBS -N MPAS_vtune_job
#PBS -A UCDV0023
#PBS -l walltime=1:00:00
#PBS -q develop
#PBS -l job_priority=regular
#PBS -j oe
#PBS -k eod
#PBS -l select=1:ncpus=64:mpiprocs=64:mem=235GB

module load vtune

export MPI_USING_VTUNE=true
export TMPDIR=/glade/derecho/scratch/$USER/$PWD
mkdir -p $TMPDIR

rm -rf vtune_output
mkdir vtune_output

# analysis runs
for analysis_type in hotspots
do
    ### Run the executable
    mpirun -n 64 vtune --collect=$analysis_type --result-dir=./vtune_output/$analysis_type --trace-mpi --return-app-exitcode -- ../../atmosphere_model
    err="${?}"
    if [ "${err}" -eq 0 ]
    then
        for x in ./vtune_output/$analysis_type*
        do
            if [ "$analysis_type" = "hotspots" ]
            then
                vtune --report=gprof-cc --result-dir=$x --format=text --report-output=./vtune_output/gprof_report_$(basename $x).txt
            elif [ "$analysis_type" = "hpc-performance" ]
            then
                vtune --report=summary -report-knob show-issues=false --result-dir=$x --format=text --report-output=./vtune_output/summary_report_$(basename $x).txt
            fi
        done
    else
        echo "failed to run $analysis_type analysis"
        break
    fi
done
