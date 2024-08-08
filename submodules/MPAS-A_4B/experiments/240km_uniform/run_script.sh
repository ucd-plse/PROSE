#!/bin/bash

set -x

costs=()
for i in 1
do
    mkdir __temp${i}

    while true; do

        if make clean > /dev/null && mpiexec -n 64 --hosts $1 ../../atmosphere_model > temp 2>&1; then
            grep -e "  1 total time" log.atmosphere.0000.out
            cat temp
            mv timing.* __temp${i}
            cd __temp${i}
            cost=$(python3 /glade/u/home/jdvanover/precimonious-w-rose/src/python/gptlparser.py | tail -n 1)
            touch COST_${cost}
            costs+=(${cost})
            cd ..
            break
        elif grep -qe "rrtmg_lw_taumol_m" temp && grep -qe "SIGSEGV" temp; then
            echo "##FAIL"
            continue
        else
            cat log.atmosphere.0000.out
            cat temp
            rm -rf __temp*
            exit 1
        fi
    done

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