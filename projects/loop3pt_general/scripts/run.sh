#!/bin/bash

run() {
    for L in {4..5}; do
	for lambda in $(seq 0.43 0.05 1.18); do
	    projects/loop3pt_general/build/release/loop3pt_general $L $lambda $1 $2 $3
	done
    done
}

run 1 1 0 || exit 1 &
run 0 1 1 || exit 1 &
run 3 1 2 || exit 1 &
run 3 2 1 || exit 1 &

wait

projects/loop3pt_general/scripts/merge.sh
