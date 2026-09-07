#!/bin/sh
# $1 = comma-separated seeds, $2 = tag
set -e
export OMP_NUM_THREADS=1
python3 sweep.py perturb --seed "$1" --tag "$2" --reps 100 --n 220
python3 sweep.py climb   --seed "$1" --tag "$2" --iters 500 --n 220 --sigma 0.03
python3 sweep.py climb   --seed "$1" --tag "${2}b" --iters 500 --n 220 --sigma 0.10 --rngseed 777
