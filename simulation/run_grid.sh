#!/bin/zsh
# Resumable: skips replicates whose JSON already exists. Four fits at a time (each fit uses four cores).
# jobs.txt has one line per replicate: "k rep" (k = 0 has 40 replicates, the other scenarios 100).
cd "$(dirname "$0")"
mkdir -p outputs
cat jobs.txt | xargs -P 4 -L 1 zsh -c 'f=$(printf "outputs/grid/rep_k%.1f_%s.json" $0 $1); [ -f "$f" ] || ./run.sh src/11_simulation_grid.py --k $0 --rep $1 >> outputs/logs_grid.txt 2>&1'
echo "DONE $(date)" >> outputs/logs_grid.txt
