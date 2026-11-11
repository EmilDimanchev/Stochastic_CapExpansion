
#=========================

Utility functions for logging and managing memory, configuring parallel workers, and running stochastic exploration with Benders and MGA

=========================#

function log_result_memory!(label::String, output::Dict)
    size_mb = round(Base.summarysize(output) / 1024^2; digits = 2)
    log_detail("Approx memory for $(label): $(size_mb) MiB")
end

function release_heavy_payload!(output::Dict)
    for key in ("SPs", "SPs_eval", "All MP outputs per iteration", "All SP outputs per iteration")
        if haskey(output, key)
            delete!(output, key)
        end
    end
    GC.gc(false)
    return nothing
end

# Every experiment entry point calls this first, so it also sets up the run's log output
# ("Log style", see src/results/Run_log.jl) before anything else prints.
function configure_parallel_workers!(settings::Dict)
    configure_run_output!(settings)
    if settings["Parallel flag"]
        desired_workers = haskey(settings, "Workers") ? settings["Workers"] : max(Sys.CPU_THREADS - 1, 1)
        wkers = nworkers() == 1 ? 0 : nworkers()
        if wkers < desired_workers
            println("Current workers: ", wkers, ". Adding ", desired_workers - wkers, " workers for parallel processing...")
            addprocs(desired_workers - wkers)
        end

        project_root = abspath(joinpath(@__DIR__, "..", ".."))
        init_expr = quote
            if !isdefined(Main, :Stochastic_CapExpansion)
                include(joinpath($project_root, "src", "Stochastic_CapExpansion.jl"))
            end
            using .Stochastic_CapExpansion
            using Revise, JuMP, Ipopt, Gurobi, HiGHS, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates, Distributions
            # One BLAS thread per worker: the workers already outnumber the cores, so a BLAS pool
            # in each one only adds threads that spin while waiting (seen as system CPU time).
            LinearAlgebra.BLAS.set_num_threads(1)
            nothing
        end
        @sync for pid in workers()
            @async remotecall_wait(Core.eval, pid, Main, init_expr)
        end

        settings["Workers"] = nworkers()
        worker_blas = remotecall_fetch(() -> LinearAlgebra.BLAS.get_num_threads(), first(workers()))
        @info("Running with parallelization using $(nworkers()) distributed workers; BLAS threads: main $(LinearAlgebra.BLAS.get_num_threads()), worker $(worker_blas); CPU threads seen: $(Sys.CPU_THREADS)")
    else
        @info("Running without parallelization")
    end
end

