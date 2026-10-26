#=========

Cut management tests: exercises add_optimality_cuts!/deactivate_cuts/reactivate_cuts
directly against a real (laptop-sized) MP, without running a full Benders loop. Checks
that a deactivate -> reactivate round trip reproduces cuts exactly (same RHS), and times
three deletion patterns at production-scale batch sizes so a regression back to any of the
slower paths would show up as a large gap: (a) constraint_by_name-per-call + scalar delete
(the original pattern), (b) cut_refs-cache-per-call + scalar delete (isolates the cache win
alone), (c) cut_refs-cache + batched delete via deactivate_cuts itself (the production
function - isolates the batch-delete win on top of the cache).

==========#

function test_cut_management(test_index = "cut_mgmt")
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)

    MP = build_planning_model(inputs, settings)

    P_s = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]
    S, F, K = length(P_s), length(P_f), length(P_k)
    G = inputs["Number of generation resources"]
    O = inputs["Number of storage resources"]
    R = G + O
    L = inputs["Number of lines"]

    x_prev = zeros(R)
    x_prev_line = zeros(L)
    coeffs = fill(1.0 / (S*F*K), S, F, K)
    SP_obj = zeros(S, F, K)
    cap_dual = zeros(R, S, F, K)
    line_dual = zeros(L, S, F, K)

    cut_refs = Dict{String, ConstraintRef}()
    n_iterations = 6
    for it in 1:n_iterations
        new_refs = add_optimality_cuts!(MP, SP_obj, cap_dual, line_dual, x_prev, x_prev_line, coeffs, inputs, settings, it, "cuttest")
        merge!(cut_refs, new_refs)
    end

    total_cuts = length(cut_refs)
    @assert total_cuts == n_iterations * S * F * K "expected $(n_iterations*S*F*K) cuts, got $(total_cuts)"
    @info("Added $(total_cuts) cuts across $(n_iterations) iterations ($(S)x$(F)x$(K) scenarios/iteration).")

    # --- Correctness: deactivate + reactivate round trip ---
    cut_archive = Dict{String, Any}()
    all_names = collect(keys(cut_refs))
    to_deactivate = all_names[1:min(50, length(all_names))]
    rhs_before = Dict(n => normalized_rhs(cut_refs[n]) for n in to_deactivate)

    deactivate_cuts(MP, to_deactivate, cut_archive, cut_refs)
    @assert length(cut_archive) == length(to_deactivate)
    @assert all(n -> !haskey(cut_refs, n), to_deactivate)
    @assert all(n -> constraint_by_name(MP, n) === nothing, to_deactivate)
    @info("Deactivated $(length(to_deactivate)) cuts; archived and removed from MP as expected.")

    reactivated = reactivate_cuts(MP, cut_archive, to_deactivate)
    @assert isempty(cut_archive)
    @assert length(reactivated) == length(to_deactivate)
    merge!(cut_refs, reactivated)

    rhs_after = Dict(n => normalized_rhs(cut_refs[n]) for n in to_deactivate)
    @assert rhs_before == rhs_after "reactivated cuts' RHS did not match the originals"
    @info("Reactivated all $(length(to_deactivate)) cuts; RHS values match the originals exactly.")

    # --- trim_cuts!: trims live + archived cuts together; both rules keep the expected cuts ---
    # Uses its own small MP so the deletions don't disturb cut_refs for the checks below.
    function make_trim_pool()
        MP_t = build_planning_model(inputs, settings)
        refs_t = Dict{String, ConstraintRef}()
        for it in 1:2
            merge!(refs_t, add_optimality_cuts!(MP_t, SP_obj, cap_dual, line_dual, x_prev, x_prev_line, coeffs, inputs, settings, it, "trimtest"))
        end
        names_t = sort(collect(keys(refs_t)))
        archive_t = Dict{String, Any}()
        deactivate_cuts(MP_t, names_t[1:20], archive_t, refs_t) # 20 archived, rest live
        # created and last_bound both ascend with position, so the rules disagree:
        # recent_binding keeps the last k, earliest keeps the first k
        stats_t = Dict{String, Tuple{Int,Int}}(n => (i, i) for (i, n) in enumerate(names_t))
        return MP_t, archive_t, stats_t, names_t
    end
    k = 50
    for (rule, pick) in (("recent_binding", names -> names[end-k+1:end]), ("earliest", names -> names[1:k]))
        MP_t, archive_t, stats_t, names_t = make_trim_pool()
        expected = Set(pick(names_t))
        n_pool = length(names_t)
        n_live_del, n_arch_del = trim_cuts!(MP_t, archive_t, stats_t, k, rule, n_pool)
        live_after = Set(name(con) for con in all_constraints(MP_t, include_variable_in_set_constraints = false) if startswith(name(con), "optimality_cut_"))
        @assert union(live_after, Set(keys(archive_t))) == expected "rule $rule kept the wrong cuts"
        @assert isempty(intersect(live_after, Set(keys(archive_t)))) "a cut ended up both live and archived"
        @assert Set(keys(stats_t)) == expected "stats not trimmed alongside the cuts"
        @assert n_live_del + n_arch_del == n_pool - k
        @assert count_live_cuts(MP_t) + length(archive_t) == k
    end
    MP_t, archive_t, stats_t, names_t = make_trim_pool()
    @assert trim_cuts!(MP_t, archive_t, stats_t, 10_000, "earliest", length(names_t)) == (0, 0) "trim should be a no-op under the cap"
    @assert count_live_cuts(MP_t) + length(archive_t) == length(names_t)
    @info("trim_cuts! trims live + archived cuts and keeps the expected ones for both rules.")

    # --- Self-healing fallback: a cut live in MP but absent from cut_refs (simulates a
    # cut carried over from a previous problem on a reused MP) should still resolve. ---
    probe_name = all_names[end]
    delete!(cut_refs, probe_name)
    con = constraint_by_name(MP, probe_name)
    @assert con !== nothing
    cut_refs[probe_name] = con # mirrors the fallback populate step in algorithm.jl
    @info("Fallback lookup for an untracked-but-live cut resolved correctly.")

    # --- Performance: constraint_by_name/scalar (a) vs cut_refs-cache/scalar (b) vs
    # cut_refs-cache/batched deactivate_cuts (c), at a production-scale batch size. n_perf
    # is intentionally much larger than the old N=300: the batch-delete win comes from
    # avoiding an O(total constraints) reindex/resync per deleted cut in HiGHS.jl/Gurobi.jl,
    # which only shows up clearly once the batch size is in the thousands (matching real
    # cuts_to_remove sizes given `Cut deactivation threshold: 8`/`Cuts retained: 50000`).
    n_perf = min(3000, length(cut_refs))
    refill_iter = Ref(n_iterations)
    refill!(needed) = begin
        iters = ceil(Int, needed / (S*F*K)) + 1
        for _ in 1:iters
            refill_iter[] += 1
            new_refs = add_optimality_cuts!(MP, SP_obj, cap_dual, line_dual, x_prev, x_prev_line, coeffs, inputs, settings, refill_iter[], "cuttest_refill")
            merge!(cut_refs, new_refs)
        end
    end

    # (a) constraint_by_name-per-call + scalar delete (the original pattern)
    perf_batch_a = collect(keys(cut_refs))[1:n_perf]
    t_a = @elapsed begin
        for n in perf_batch_a
            con = constraint_by_name(MP, n)
            delete(MP, con)
        end
    end
    @info("(a) constraint_by_name + scalar delete: removed $(n_perf) cuts in $(round(t_a; digits=4))s")
    for n in perf_batch_a
        delete!(cut_refs, n)
    end

    # (b) cut_refs cache + scalar delete (isolates the cache win alone, no batching)
    refill!(n_perf)
    perf_batch_b = collect(keys(cut_refs))[1:n_perf]
    t_b = @elapsed begin
        for n in perf_batch_b
            delete(MP, cut_refs[n])
        end
    end
    @info("(b) cut_refs cache + scalar delete: removed $(n_perf) cuts in $(round(t_b; digits=4))s")
    for n in perf_batch_b
        delete!(cut_refs, n)
    end

    # (c) cut_refs cache + batched delete, via the production deactivate_cuts function
    refill!(n_perf)
    perf_batch_c = collect(keys(cut_refs))[1:n_perf]
    cut_archive_c = Dict{String, Any}()
    t_c = @elapsed deactivate_cuts(MP, perf_batch_c, cut_archive_c, cut_refs)
    @info("(c) cut_refs cache + batched delete (deactivate_cuts): removed $(n_perf) cuts in $(round(t_c; digits=4))s")

    speedup_cache = t_a / max(t_b, 1e-9)
    speedup_batch = t_b / max(t_c, 1e-9)
    speedup_total = t_a / max(t_c, 1e-9)
    @info("[$(test_index)] Speedup from cut_refs cache alone: $(round(speedup_cache; digits=1))x; from batching on top of the cache: $(round(speedup_batch; digits=1))x; combined: $(round(speedup_total; digits=1))x")

    return (n_perf = n_perf, scalar_by_name_time = t_a, scalar_cached_time = t_b, batched_time = t_c,
            speedup_cache = speedup_cache, speedup_batch = speedup_batch, speedup_total = speedup_total,
            total_cuts = total_cuts)

end

#=========

Cut trim rule comparison: runs the same laptop risk-pareto + MGA exploration (same seed, so the
same MGA vectors and budget schedule) once with archive trimming effectively off (baseline), then
once per trim rule ("recent_binding", "earliest") with a small "Cuts retained". Objectives/gaps
should agree across configs; the signal for which rule keeps the more useful cuts is the number
of Benders iterations needed (a worse cut set needs more). Solve *times* are reported but are not
meaningful at laptop scale (see CLAUDE.md traps), so don't read a performance win from them.

==========#

function test_cut_trim_rules(test_index = "cut_trim"; cuts_retained_small::Int = 20000, n_vectors::Int = 3)
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    base_settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, base_settings)
    configure_parallel_workers!(base_settings)
    SPs = build_all_subproblems(inputs, base_settings)

    configs = [("baseline", 10^9, "recent_binding"),
               ("recent_binding", cuts_retained_small, "recent_binding"),
               ("earliest", cuts_retained_small, "earliest")]
    summaries = Dict{String, Any}()
    for (label, n_keep, rule) in configs
        settings = deepcopy(base_settings)
        settings["Cuts retained"] = n_keep
        settings["Iterations"] = n_vectors # more than one trim event
        settings["Cut trim rule"] = rule
        settings["Cut deactivation strategy"] = "in benders"
        results_folder = joinpath("outputs", "Test_"*string(test_index), label)
        summary_folder = joinpath(results_folder, "Summary")
        mkpath(summary_folder)
        run_log = []
        t = @elapsed run_stochastic_exploration_risk_pareto(SPs, inputs, settings, joinpath(results_folder, "Pareto"), summary_folder;
            budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
            mapping_settings = Dict("Mapping" => false),
            run_settings = Dict("Summary name" => "pareto", "Run log" => run_log))
        mga_runs = filter(r -> startswith(r.label, "Random_"), run_log)
        summaries[label] = (total_iterations = sum(r.iterations for r in run_log), mga_iterations = sum(r.iterations for r in mga_runs), mp_time = sum(r.mp_time for r in run_log), wall_time = t, max_archive_at_run_start = maximum(r.archive_before for r in run_log), run_log = run_log)
        # Invariant: trimming happens at vector/risk-weight starts, so the archive going *into*
        # the first run of each vector must already be within the cap.
        if n_keep < 10^9
            first_runs = filter(r -> occursin(r"^Random_\d+(_Budget_1)?$", r.label) || startswith(r.label, "Risk_Weight_"), run_log)
            @assert all(r -> r.archive_before <= n_keep, first_runs) "live + archived cuts exceeded Cuts retained at a vector boundary for $label"
        end
    end

    for (label, s) in summaries
        @info("[$(test_index)] $label: total Benders iterations = $(s.total_iterations), MGA iterations = $(s.mga_iterations), MP time = $(round(s.mp_time; digits=3))s, wall = $(round(s.wall_time; digits=1))s, max archive at run start = $(s.max_archive_at_run_start)")
    end
    return summaries
end
