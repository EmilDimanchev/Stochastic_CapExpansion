#!/bin/bash
# GridFill run on Della with the cut archive export on. Writes
#   outputs/Test_<index>/Summary/cut_archive_GridFill.h5
# next to Summary_GridFill.csv. Submit from the repo root: sbatch experiments/cut_export/export_cuts_della.sh
#
# Before the first run:
#   - On a login node (compute nodes have no internet), install the new HDF5 dependency:
#       module add julia/1.9.1 && julia --project=. -e 'using Pkg; Pkg.instantiate()'
#   - Pull this branch so inputs/Inputs_2z_7tech/settings.yaml has the "Cut export settings" block.
#     A missing block stops the run at startup, not at the end.
#
# What the export adds to the job. These are scaled from the laptop run (5 scenarios, 294 raw cuts
# per scenario: pruning 3.7 s), assuming 125 scenarios and about 5,300 raw cuts per scenario. Check
# them against the "Cut export: pruning and writing" log line of the first real run.
#   - Time:   about 30 min per budget step (one step unless the file is over 25 MB), and up to about
#             30 min of validation, about the same work as the sampler's own per-batch pruning.
#   - Memory: under 1 GB on the main process: cuts about 75 MB, held-out blends about 90 MB, plus
#             temporaries. The 2 GB margin below is for safety.
#   - Disk:   a raw copy of about 100 MB while the export runs; removed on success, kept on failure.
#             To prune it again without the model:
#               export_cut_archive_from_raw(raw_path, summary_folder, settings["Cut export settings"]; file_stem = "cut_archive_GridFill")
#
#SBATCH --job-name=GridFill_cut_export
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2100M                 # 2G per CPU as in the Pareto runs, plus about 2 GB in total for the export
#SBATCH --time=<run time + 1:30:00>         # set this: your usual GridFill time plus the export estimate above
#SBATCH --output="cut_export_%j.out"
#SBATCH --error="cut_export_%j.err"

module add julia/1.9.1
module add gurobi/12.0.0

julia --project=. -e '
    include("experiments/MGA_Tests.jl")
    gridfill_test_della("cut_export_1"; cut_export_overrides = Dict("Export cuts" => true))
'
