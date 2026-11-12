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
# What the export adds to the job. The first 2-zone run (125 scenarios × 5 blocks, about 3,200 raw
# cuts per subproblem, 2,661 plans) took 2.9 h to prune and 28 min to validate, all on the main
# process. Pruning and the validation replay now run one subproblem per remote call on the workers,
# so expect minutes with 64 workers; check the "pruned in ... s" and "repair rounds" columns of the
# "Cut export budget step" line and the "Cut export: pruning and writing" line.
#   - Time:   repair rounds dominate: each evaluates every cut of a subproblem at "Held-out blends"
#             points, up to "Max repair rounds" times (the first run averaged about 30 of 40).
#             A second budget step (file over 25 MB) repeats the pruning.
#   - Memory: main holds the cuts (about 240 MB for 2M cuts); each worker rebuilds the blend sets
#             (about 95 MB) and a few hundred KB of cuts per call. Cleared when the export ends.
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
