#=========

Cut pool tests (src/model/benders/cut_pool.jl), run against a real laptop-sized MP without a
full Benders loop:

  test_cut_management - archive -> restore round trip reproduces cuts exactly, trim_cuts!
                        keeps the expected cuts under both rules, and times three deletion
                        patterns at production batch sizes: (a) constraint_by_name + scalar
                        delete (the original pattern), (b) cached refs + scalar delete, (c)
                        the pool's _archive! (archive copy + batched delete).
  test_cut_pool       - the vectorised archive check matches JuMP's own evaluation (before
                        and after compaction), most_violated orders and caps correctly,
                        separate! leaves nothing violated, and seed!/prune! manage the
                        per-vector binding sets.

==========#

# Laptop MP with risk terms plus `n_iterations` rounds of random (but realistic-signed) cuts,
# registered in a fresh pool. Shared by the cut tests below.
function _cut_test_setup(n_iterations::Int; case = "cuttest", seed = 1)
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    settings = load_settings(inputs_folder)
    settings["Risk aversion flag"] = true # so cvar_tail_cuts_* (and u, ζ) are exercised too
    inputs = load_input_data(inputs_folder, settings)
    MP = build_planning_model(inputs, settings)
    pool = CutPool(deepcopy(settings))

    rng = MersenneTwister(seed)
    for it in 1:n_iterations
        add_cuts!(pool, _add_random_cuts!(MP, inputs, settings, it, case, rng))
    end
    S, F, K = size(inputs["Demand scenario probabilities"], 1), size(inputs["Gas price scenario probabilities"], 1), size(inputs["Weather scenario probabilities"], 1)
    return MP, pool, inputs, settings, S*F*K
end

# One Benders iteration's worth of cuts (both kinds) from random SP outputs; returns their refs.
function _add_random_cuts!(MP, inputs, settings, it, case, rng)
    S, F, K = size(inputs["Demand scenario probabilities"], 1), size(inputs["Gas price scenario probabilities"], 1), size(inputs["Weather scenario probabilities"], 1)
    R = inputs["Number of generation resources"] + inputs["Number of storage resources"]
    L = inputs["Number of lines"]
    sp_obj = settings["Scaling factor cost"] .* (1 .+ rand(rng, S, F, K))
    cap_dual = -settings["Scaling factor cost"] .* rand(rng, R, S, F, K) ./ 1e3 # op cost falls with capacity
    line_dual = -settings["Scaling factor cost"] .* rand(rng, L, S, F, K) ./ 1e3
    x_prev = rand(rng, R) .* inputs["Capacity upper bounds"]
    return add_optimality_cuts!(MP, sp_obj, cap_dual, line_dual, x_prev, rand(rng, L), ones(S, F, K), inputs, settings, it, case)
end

_cut_signature(con) = (obj = constraint_object(con); (Dict(obj.func.terms), MOI.constant(obj.set) - obj.func.constant))

function test_cut_management(test_index = "cut_mgmt")
    n_iterations = 6
    MP, pool, inputs, settings, n_scen = _cut_test_setup(n_iterations)
    total_cuts = length(pool.live)
    @assert total_cuts == 2 * n_iterations * n_scen "expected $(2*n_iterations*n_scen) cuts, got $(total_cuts)"
    @assert count_live_cuts(MP) == total_cuts
    @info("Added $(total_cuts) cuts across $(n_iterations) iterations ($(n_scen) scenarios/iteration, optimality + CVaR tail).")

    # --- Correctness: archive + restore round trip ---
    to_archive = sort(collect(keys(pool.live)))[1:50]
    before = Dict(n => _cut_signature(pool.live[n]) for n in to_archive)
    @assert _archive!(pool, MP, to_archive) == length(to_archive)
    @assert n_archived(pool) == length(to_archive) && all(n -> !haskey(pool.live, n), to_archive)
    @assert all(n -> constraint_by_name(MP, n) === nothing, to_archive)
    _restore!(pool, MP, [pool.colof[n] for n in to_archive])
    @assert n_archived(pool) == 0 && count_live_cuts(MP) == total_cuts
    @assert all(n -> _cut_signature(pool.live[n]) == before[n], to_archive) "restored cuts differ from the originals"
    @info("Archived and restored $(length(to_archive)) cuts; coefficients and RHS match exactly.")

    # --- trim_cuts!: trims live + archived cuts together; both rules keep the expected cuts ---
    k = 50
    function make_trim_pool(rule)
        MP_t, pool_t, _, _, _ = _cut_test_setup(1; case = "trimtest")
        names_t = sort(collect(keys(pool_t.live)))
        _archive!(pool_t, MP_t, names_t[1:20]) # 20 archived, rest live
        # created and last_bound both ascend with position, so the rules disagree:
        # recent_binding keeps the last k, earliest keeps the first k
        pool_t.stats = Dict(n => (i, i) for (i, n) in enumerate(names_t))
        pool_t.clock = length(names_t)
        pool_t.cfg["Cuts retained"], pool_t.cfg["Trim rule"] = k, rule
        return MP_t, pool_t, names_t
    end
    for (rule, pick) in (("recent_binding", names -> names[end-k+1:end]), ("earliest", names -> names[1:k]))
        MP_t, pool_t, names_t = make_trim_pool(rule)
        n_live_del, n_arch_del = trim_cuts!(MP_t, pool_t)
        live_after = Set(keys(pool_t.live))
        @assert union(live_after, Set(keys(pool_t.colof))) == Set(pick(names_t)) "rule $rule kept the wrong cuts"
        @assert isempty(intersect(live_after, Set(keys(pool_t.colof)))) "a cut ended up both live and archived"
        @assert Set(keys(pool_t.stats)) == Set(pick(names_t)) "stats not trimmed alongside the cuts"
        @assert n_live_del + n_arch_del == length(names_t) - k
        @assert count_live_cuts(MP_t) == length(live_after) && length(pool_t.names) == n_archived(pool_t) "trim should compact the archive"
    end
    MP_t, pool_t, names_t = make_trim_pool("earliest")
    pool_t.cfg["Cuts retained"] = 10_000
    @assert trim_cuts!(MP_t, pool_t) == (0, 0) "trim should be a no-op under the cap"
    @info("trim_cuts! trims live + archived cuts and keeps the expected ones for both rules.")

    # --- Performance: constraint_by_name/scalar (a) vs cached refs/scalar (b) vs the pool's
    # batched _archive! (c), at a production-scale batch size. The batch-delete win comes
    # from avoiding an O(total constraints) reindex/resync per deleted cut in HiGHS.jl/
    # Gurobi.jl, which only shows up clearly once the batch size is in the thousands.
    n_perf = 3000
    refill_iter = Ref(n_iterations)
    rng = MersenneTwister(2)
    fresh_batch!() = begin
        refs = Dict{String, ConstraintRef}()
        while length(refs) < n_perf
            refill_iter[] += 1
            merge!(refs, _add_random_cuts!(MP, inputs, settings, refill_iter[], "cuttest_refill", rng))
        end
        return collect(keys(refs))[1:n_perf], refs
    end

    names_a, _ = fresh_batch!()
    t_a = @elapsed for n in names_a
        delete(MP, constraint_by_name(MP, n))
    end
    @info("(a) constraint_by_name + scalar delete: removed $(n_perf) cuts in $(round(t_a; digits=4))s")

    names_b, refs_b = fresh_batch!()
    t_b = @elapsed for n in names_b
        delete(MP, refs_b[n])
    end
    @info("(b) cached refs + scalar delete: removed $(n_perf) cuts in $(round(t_b; digits=4))s")

    names_c, refs_c = fresh_batch!()
    add_cuts!(pool, refs_c)
    t_c = @elapsed _archive!(pool, MP, names_c)
    @info("(c) _archive! (archive copy + batched delete): removed $(n_perf) cuts in $(round(t_c; digits=4))s")

    speedup_cache = t_a / max(t_b, 1e-9)
    speedup_batch = t_b / max(t_c, 1e-9)
    speedup_total = t_a / max(t_c, 1e-9)
    @info("[$(test_index)] Speedup from cached refs alone: $(round(speedup_cache; digits=1))x; from batching on top: $(round(speedup_batch; digits=1))x; combined: $(round(speedup_total; digits=1))x")

    return (n_perf = n_perf, scalar_by_name_time = t_a, scalar_cached_time = t_b, batched_time = t_c,
            speedup_cache = speedup_cache, speedup_batch = speedup_batch, speedup_total = speedup_total,
            total_cuts = total_cuts)
end

function test_cut_pool(test_index = "cut_pool")
    MP, pool, inputs, settings, n_scen = _cut_test_setup(4)
    all_names = sort(collect(keys(pool.live)))
    sigs = Dict(n => _cut_signature(pool.live[n]) for n in all_names)

    # Manual violation of cut n at variable values val (VariableRef -> Float64)
    manual_viol(n, val) = ((terms, rhs) = sigs[n]; (rhs - sum(c * val[v] for (v, c) in terms)) / max(1.0, abs(rhs)))
    function check_violations(label)
        v = randn(MersenneTwister(7), length(pool.vars)) .* 100
        val = Dict(zip(pool.vars, v))
        viol = archive_violations(pool, v)
        for (n, j) in pool.colof
            @assert isapprox(viol[j], manual_viol(n, val); rtol = 1e-10, atol = 1e-10) "vectorised violation of $n differs ($label)"
        end
        @assert all(isinf, viol[.!pool.alive]) "dead columns must not count as violated ($label)"
        return v, viol
    end

    # --- Vectorised violation check matches a per-cut evaluation, with dead columns present
    # (after a partial restore) and after compaction ---
    _archive!(pool, MP, all_names)
    @assert isempty(pool.live) && n_archived(pool) == length(all_names)
    check_violations("fresh archive")
    _restore!(pool, MP, 1:length(all_names) ÷ 3) # a third restored: dead columns, no compaction yet
    @assert length(pool.names) == length(all_names) "a third dead should not trigger compaction"
    check_violations("with dead columns")
    _restore!(pool, MP, [pool.colof[n] for n in collect(keys(pool.colof))[1:length(all_names) ÷ 4]]) # now dead > half
    @assert length(pool.names) == n_archived(pool) "dead columns past half should trigger compaction"
    v, viol = check_violations("after compaction")
    @info("Vectorised archive check matches per-cut evaluation for $(n_archived(pool)) cuts (with dead columns and after compaction).")

    # --- most_violated: above tol, most violated first, capped ---
    tol, cap = 1e-6, 25
    picked = most_violated(pool, v, tol, cap)
    n_above = count(>(tol), viol)
    @assert length(picked) == min(cap, n_above)
    @assert issorted(viol[picked]; rev = true) "most_violated is not most-violated-first"
    @assert isempty(picked) || minimum(viol[picked]) >= maximum(viol[setdiff(findall(>(tol), viol), picked)]; init = -Inf) "most_violated skipped a more violated cut"
    @info("most_violated picked $(length(picked)) of $(n_above) violated cuts, in order.")

    # --- separate!: after it returns, no archived cut is violated (unless the round limit hit) ---
    pool.cfg["Max rounds"] = 50
    output, sep = separate!(pool, MP, settings, settings["Risk aversion weight"])
    @assert sep.rounds < pool.cfg["Max rounds"] "separation should finish well inside 50 rounds here"
    @assert maximum(archive_violations(pool, value.(pool.vars)); init = -Inf) <= pool.cfg["Separation tolerance"] "separate! left an archived cut violated"
    @assert all(s -> s.restored <= pool.cfg["Cuts per round"] * s.rounds, [sep])
    @info("separate! restored $(sep.restored) cuts in $(sep.rounds) rounds; nothing archived is violated at the final MP solution.")

    # --- prune!: records binding cuts in bound_cur, archives non-binding ones after `window` ---
    window = 2
    prune!(pool, MP, 1, window) # every live cut is new to this solve: grace only
    binding = Set(n for (n, con) in pool.live if dual(con) >= pool.cfg["Binding tolerance"])
    @assert binding ⊆ pool.bound_cur
    @assert !isempty(binding) "expected some cut to bind"
    prune!(pool, MP, 1 + window, window)
    @assert all(n -> n in binding || !haskey(pool.live, n), keys(sigs)) "a cut that never bound survived past the window"
    @assert binding ⊆ Set(keys(pool.live)) "a binding cut was archived"
    @info("prune! kept the $(length(binding)) binding cuts and archived the rest after $window iterations.")

    # --- seed!: live set = this vector's binding set within a vector, previous vector's on a new one ---
    pool.group = 1
    pool.bound_cur = Set(all_names[1:30])
    seed!(pool, MP, 1)
    @assert Set(keys(pool.live)) == Set(all_names[1:30]) "same-vector seed should be the vector's binding set"
    @assert isempty(pool.last_bound_iter) "seed! must reset per-solve grace"
    seed!(pool, MP, 2)
    @assert pool.bound_prev == Set(all_names[1:30]) && isempty(pool.bound_cur) "binding sets did not roll over on a new vector"
    @assert Set(keys(pool.live)) == pool.bound_prev "new-vector seed should be the previous vector's binding set"
    @assert count_live_cuts(MP) == length(pool.live) && length(pool.live) + n_archived(pool) == length(all_names) "a cut was lost or duplicated"
    @info("[$(test_index)] seed! picks the right live set within and across vectors; all cut pool checks passed.")
    return true
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
        settings["Cut settings"]["Cuts retained"] = n_keep
        settings["Iterations"] = n_vectors # more than one trim event
        settings["Cut settings"]["Trim rule"] = rule
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


#=========

Cut pool A/B: runs the laptop risk-pareto + MGA exploration (anchors, then `n_vectors` MGA
vectors with budget tightening and wing caps) on this checkout, writing the Summary CSV and the
per-solve run log to outputs/Test_<test_index>/new. Pass `old_folder` - the same run made with
the pre-pool code (run_cut_pool_ab_case copied into a checkout of the old commit; its run log
only has label/iterations/mp_time/archive_before columns) - to compare per solve:

  - final capacity mix and Transformed / CVaR / EV / Inv values must agree within the 0.1%
    convergence tolerance
  - SP rounds, MP solves, live cuts at the first MP solve, MP time and separation time are
    reported; solves whose SP rounds rise by more than 3, or whose separation time exceeds
    10% of their MP time, are flagged

Laptop solve times are overhead-bound (see CLAUDE.md traps); read iterations / SP rounds / live
cut counts, not seconds. Writes outputs/Test_<test_index>/comparison.csv.

==========#

function run_cut_pool_ab_case(folder::String; n_vectors::Int = 3, settings_overrides::Dict = Dict())
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    settings = load_settings(inputs_folder)
    merge!(settings, settings_overrides)
    settings["Iterations"] = n_vectors
    inputs = load_input_data(inputs_folder, settings)
    configure_parallel_workers!(settings)
    SPs = build_all_subproblems(inputs, settings)
    summary_folder = joinpath(folder, "Summary")
    mkpath(summary_folder)
    run_log = []
    run_stochastic_exploration_risk_pareto(SPs, inputs, settings, joinpath(folder, "Pareto"), summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
        mapping_settings = Dict("Mapping" => false),
        run_settings = Dict("Summary name" => "pareto", "Run log" => run_log))
    CSV.write(joinpath(folder, "run_log.csv"), DataFrame(run_log))
    return run_log
end

function test_cut_pool_ab(test_index = "cut_pool_ab"; n_vectors::Int = 3, old_folder::Union{String, Nothing} = nothing, run_new::Bool = true)
    root = joinpath("outputs", "Test_"*string(test_index))
    new_folder = joinpath(root, "new")
    run_new && run_cut_pool_ab_case(new_folder; n_vectors = n_vectors)
    old_folder === nothing && return nothing
    comparison = compare_cut_pool_ab(old_folder, new_folder)
    CSV.write(joinpath(root, "comparison.csv"), comparison)
    return comparison
end

function compare_cut_pool_ab(old_folder::String, new_folder::String; tol = 1e-3)
    read_summary(f) = CSV.read(joinpath(f, "Summary", "Summary_pareto.csv"), DataFrame)
    sum_old, sum_new = read_summary(old_folder), read_summary(new_folder)
    log_old = CSV.read(joinpath(old_folder, "run_log.csv"), DataFrame)
    log_new = CSV.read(joinpath(new_folder, "run_log.csv"), DataFrame)
    cap_cols = names(sum_old)[1:findfirst(==("Index"), names(sum_old))-1]

    # Transformed cost, with each run's own anchor-based transform (as in the driver)
    function transformed(df)
        row(w) = df[findfirst(==("Risk_Weight_"*string(w)), df.Index), :]
        a0, a1 = row(0.0), row(1.0)
        k_cvar = 1 / (a1.CVaR_OpCost - a0.CVaR_OpCost)
        k_sys = 1 / (a0.EV_SystemCost - a1.EV_SystemCost)
        return k_cvar .* df.CVaR_OpCost .+ k_sys .* df.EV_SystemCost
    end
    sum_old.Transformed, sum_new.Transformed = transformed(sum_old), transformed(sum_new)

    reldiff(a, b) = abs(b - a) / max(abs(a), 1.0)
    rows = []
    for r_new in eachrow(sum_new)
        i_old = findfirst(==(r_new.Index), sum_old.Index)
        i_old === nothing && (@warn("Solve $(r_new.Index) has no old counterpart"); continue)
        r_old = sum_old[i_old, :]
        l_old = log_old[findfirst(==(r_new.Index), log_old.label), :]
        l_new = log_new[findfirst(==(r_new.Index), log_new.label), :]
        diffs = (cap = maximum(reldiff(r_old[c], r_new[c]) for c in cap_cols),
                 transformed = reldiff(r_old.Transformed, r_new.Transformed),
                 cvar = reldiff(r_old.CVaR_OpCost, r_new.CVaR_OpCost),
                 ev = reldiff(r_old.EV_OpCost, r_new.EV_OpCost),
                 inv = reldiff(r_old.Investment_Cost, r_new.Investment_Cost))
        flags = String[]
        maximum(values(diffs)) > tol && push!(flags, "result_mismatch")
        l_new.sp_rounds - l_old.iterations > 3 && push!(flags, "sp_rounds_up")
        l_new.sep_time > 0.1 * l_new.mp_time && push!(flags, "sep_time_gt_10pct")
        push!(rows, (solve = r_new.Index,
                     sp_rounds_old = l_old.iterations, sp_rounds_new = l_new.sp_rounds, mp_solves_new = l_new.mp_solves,
                     live_at_start_old = l_old.archive_before, live_at_start_new = l_new.live_at_start,
                     mp_time_old = l_old.mp_time, mp_time_new = l_new.mp_time, sep_time_new = l_new.sep_time,
                     maxreldiff_cap = diffs.cap, reldiff_transformed = diffs.transformed, reldiff_cvar = diffs.cvar,
                     reldiff_ev = diffs.ev, reldiff_inv = diffs.inv, flags = join(flags, ";")))
    end
    df = DataFrame(rows)
    n_flag = count(!isempty, df.flags)
    @info("Cut pool A/B: $(nrow(df)) solves compared, $n_flag flagged; SP rounds old/new = $(sum(df.sp_rounds_old))/$(sum(df.sp_rounds_new)), live cuts at first MP old/new = $(sum(df.live_at_start_old))/$(sum(df.live_at_start_new))")
    return df
end
