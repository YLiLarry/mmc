#!/usr/bin/env bash

for d in 32 64; do
    for b in `seq 21 23`; do
        # for e in `seq 15 18`; do
            e=$(($b - 4))
            if [[ b -gt e ]]; then
                echo -d ${d} -b ${b} -e ${e}
                ./eugene-eric-larry-research-do-not-kill -d ${d} -b ${b} -e ${e} > log/d${d}_b${b}_e${e}.txt 2> log/d${d}_b${b}_e${e}_stderr.txt
            fi
        # done
    done
done
