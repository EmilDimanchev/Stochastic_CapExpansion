#!/bin/bash
#SBATCH --job-name=Stability_test
#SBATCH --nodes=1                           # node count
#SBATCH --ntasks=1                          # total number of tasks across all nodes
#SBATCH --cpus-per-task=32                   # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=1G                    # memory per cpu-core
#SBATCH --time=2:00:00                     # total run time limit (HH:MM:SS)
#SBATCH --output="test.out"
#SBATCH --mail-type=end                    # notifications for job done & fail
#SBATCH --mail-user=ml6802@princeton.edu  # send-to address


module add julia/1.9.1
module add gurobi/12.0.0
julia --project="/home/ml6802/Stochastic_CapExpansion" -p64 stability_test.jl

date