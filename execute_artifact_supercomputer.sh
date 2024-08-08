#!/bin/bash

PBS_ACCOUNT=$1

# set environment
source $(git rev-parse --show-toplevel)/scripts/set_env.sh

# set cwd to repo root
cd ${PROSE_REPO_PATH}

set -x

# compile plugin
pushd src/cpp
make
popd

# set paths in setup.ini
for x in ADCIRC_4B/ MOM6_4B/ MPAS-A_4B/ MPAS-A_4C/
do
    sed -i "s|<REPLACE_WITH_PATH_TO_REPO_ROOT>|$(realpath submodules/$x)|g" $(find submodules/$x -name setup.ini)
done

# unzip MPAS-A input files
for x in MPAS-A_4B/ MPAS-A_4C/
do
    pushd submodules/$x/experiments/240km_uniform
    unzip x1.10242.init.nc.zip
    popd
done

# submit experiment jobs
for x in submodules/ADCIRC_4B/experiments/beaufort/ submodules/MOM6_4B/experiments/benchmark_zstar/ submodules/MPAS-A_4B/experiments/240km_uniform/ submodules/MPAS-A_4C/experiments/240km_uniform/
do
    pushd $x
    qsub -A ${PBS_ACCOUNT} -l walltime=12:00:00 -v PROSE_REPO_PATH=${PROSE_REPO_PATH} search_job.sh
    popd
done

# wait for minimal amount of time to pass
sleep 12h

# check every 5 minutes for end of all running experiments
while qstat -u ${USER} > /dev/null 2>&1
do
    sleep 5m
done

# generate figures
for x in submodules/ADCIRC_4B/experiments/beaufort/ submodules/MOM6_4B/experiments/benchmark_zstar/ submodules/MPAS-A_4B/experiments/240km_uniform/ submodules/MPAS-A_4C/experiments/240km_uniform/
do
    pushd $x
    prose_check_progress.py > prose_logs/search_log.txt
    python3 generate_figures.py > /dev/null 2>&1
    mv *.html ${PROSE_REPO_PATH}
    popd
done