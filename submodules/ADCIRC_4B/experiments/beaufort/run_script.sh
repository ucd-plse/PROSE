#!/bin/bash

set -x

costs=()
for i in 1
do
    mkdir __temp${i}
    if mpiexec -n 128 --hosts $1 ../../work/padcirc; then
        mv timing.* __temp${i}
        cd __temp${i}
        cost=$(python3 /glade/u/home/jdvanover/precimonious-w-rose/src/python/gptlparser.py | tail -n 1)
        touch COST_${cost}
        costs+=(${cost})
        cd ..
    else
        rm __temp${i} -rf
        exit 1
    fi
done

# Sort the array
sorted_costs=($(for i in "${costs[@]}"; do echo "$i"; done | sort -n))

# Find the index of the median value
array_length=${#sorted_costs[@]}
median_index=$((array_length / 2))
median_cost=${sorted_costs[${median_index}]}

mv $(dirname $(find -name "COST_${median_cost}"))/timing.* .
rm __temp*  -rf
exit 0