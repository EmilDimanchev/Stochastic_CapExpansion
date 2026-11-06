import HDF5
using LinearAlgebra, Printf, Random, Serialization, SparseArrays, Statistics

#==========================================================================================

Pruned Benders cut archive for the Grid Options Explorer

The dashboard builds blends x̄ = Σ w_i x_i of solved plans and needs lower bounds on each
scenario's operating cost at those blends:
    LB_s(x̄) = max_k (α_k + β_k·x̄),  over the cuts k of scenario s.
This file picks the few cuts that matter for such blends and writes them to one HDF5 file.

Flow (export_cut_archive):
  1. collect   pool cuts + final-solve cuts, unscaled and grouped by scenario
  2. save      a temporary raw copy, so a failure here can't lose the cuts of a long run
  3. prune     one scenario at a time:
                 dedupe -> keep cuts active at probe blends -> check held-out blends -> repair
  4. write     the HDF5 file, shrinking the probe set if it is over budget
  5. validate  reproduce the logged GridFill bounds, check validity and tightness at the plans

Units: α in $, β in $/MW, x in MW. θ_s is the annual, undiscounted operating cost of scenario
s (variable generation + non-served energy + capacity-reserve slack; no investment).

File layout, as h5py shows it (K cuts, n first-stage variables, N plans). Julia stores arrays
column-major, so a Julia n×K matrix reads as K×n in h5py.
  cuts/scenario    K, int16, 1-based, index into meta/scenario_label and meta/scenario_probability
  cuts/alpha       K, float64
  cuts/beta        K×n, float32 or float64 (or cuts/beta_csr/{indptr,indices,data}, 0-based,
                   when sparse; meta/beta_format says which)
  cuts/iteration   K, int16, Benders iteration that generated the cut (0 for a final solve)
  cuts/source      K, int16, 1-based, index into meta/sources: the solve that generated the cut
  cuts/batch       K, int16, GridFill batch that generated the cut (0 for anything before GridFill)
  vars/{name,resource,zone,units}   n strings, column order of beta and x
  plans/id         N strings, the Index column of the Summary CSV
  plans/x          N×n, float32 or float64
  meta/...         run id, probabilities, CVaR α, units, pruning settings, counts, gap statistics

==========================================================================================#

const CUT_EXPORT_CHUNK_PROBES = 500   # probes evaluated at once; bounds the K × chunk temporary

#---- 1. Collect ----------------------------------------------------------------------------

# Cuts per scenario. Within a scenario the order is: pool cuts, then final solves in bank order.
struct CutTable
    alpha::Vector{Vector{Float64}}   # per scenario: K_s intercepts, $
    beta::Vector{Matrix{Float64}}    # per scenario: n × K_s slopes, $/MW
    source::Vector{Vector{Int}}      # per scenario: index into CutExportInput.sources
    iteration::Vector{Vector{Int}}   # per scenario: Benders iteration (0 for a final solve)
    batch::Vector{Vector{Int}}       # per scenario: GridFill batch (0 for anything before GridFill)
end

n_table_cuts(t::CutTable) = sum(length, t.alpha)

# Everything the export needs, copied out of the run. Only plain data, so it can be saved
# with Serialization and pruned again without the model.
struct CutExportInput
    run_id::String
    cuts::CutTable
    sources::Vector{String}
    skipped_rows::Int                # pool optimality rows with no alpha term (non-optimal subproblem)
    plan_ids::Vector{String}
    plan_x::Matrix{Float64}          # n × N
    plan_θ::Matrix{Float64}          # S × N, from subproblems solved at plan_x
    plan_ok::Vector{Bool}            # false if any of the plan's subproblems failed
    plan_batch::Vector{Int}          # GridFill batch of the plan (0 for anything before GridFill)
    var_name::Vector{String}
    var_resource::Vector{String}
    var_zone::Vector{String}
    var_units::Vector{String}
    p::Vector{Float64}               # S scenario probabilities
    scenario_label::Vector{String}   # D_<s>_F_<f>_W_<k>, as in the Summary's OperatingCost_ columns
    cvar_alpha::Float64
    seed::Int
    scaling_factor_cost::Float64
    picks::Vector{GridFillPick}      # empty if GridFill didn't run
end

# Cut names are optimality_cut_<case>_<iteration>_[s,f,k]; see add_optimality_cuts!.
function parse_cut_name(name::AbstractString)
    m = match(r"^optimality_cut_(.*)_(\d+)_\[", name)
    m === nothing && error("Unexpected cut name \"$name\"")
    return String(m.captures[1]), parse(Int, m.captures[2])
end

# Group pool cuts and final-solve cuts by scenario. `pool_cuts` comes from pool_optimality_cuts.
function collect_cut_table(pool_cuts, bank::AbstractVector, plan_batch::AbstractVector{Int})
    S, n = length(bank[1].θ), length(bank[1].x)
    alpha = [Float64[] for _ in 1:S]
    slopes = [Vector{Float64}[] for _ in 1:S]
    source, iteration, batch = ([Int[] for _ in 1:S] for _ in 1:3)
    sources = String[]
    source_id = Dict{String, Int}()
    id_of(label) = get!(source_id, label) do
        push!(sources, label)
        length(sources)
    end
    function add!(s, a, g, label, iter, b)
        push!(alpha[s], a); push!(slopes[s], g); push!(source[s], id_of(label))
        push!(iteration[s], iter); push!(batch[s], b)
    end
    for cut in pool_cuts
        case, iter = parse_cut_name(cut.name)
        add!(cut.scenario, cut.a, cut.g, case, iter, 0)
    end
    # The same tangent add_final_solve! builds: θ_s(x) >= θ_s(x_j) + g'(x - x_j)
    for (fs, b) in zip(bank, plan_batch)
        (fs.ok && has_gradients(fs)) || continue
        for s in 1:S
            g = fs.β[:, s]
            add!(s, fs.θ[s] - dot(g, fs.x), g, fs.name, 0, b)
        end
    end
    beta = [isempty(slopes[s]) ? zeros(n, 0) : reduce(hcat, slopes[s]) for s in 1:S]
    return CutTable(alpha, beta, source, iteration, batch), sources
end

# Input from toy or real parts; the inputs/settings wrapper is below.
function build_cut_export_input(pool_cuts, skipped_rows::Int, bank::AbstractVector, picks::AbstractVector;
                                run_id::String, var_name, var_resource, var_zone, var_units, scenario_label,
                                cvar_alpha::Real, seed::Integer, scaling_factor_cost::Real)
    isempty(bank) && error("Cut export: no final solves were recorded.")
    batch_of = Dict(pk.name => pk.batch for pk in picks)
    plan_batch = [get(batch_of, fs.name, 0) for fs in bank]
    table, sources = collect_cut_table(pool_cuts, bank, plan_batch)
    return CutExportInput(run_id, table, sources, skipped_rows,
                          [fs.name for fs in bank], reduce(hcat, [fs.x for fs in bank]), reduce(hcat, [fs.θ for fs in bank]),
                          [fs.ok for fs in bank], plan_batch, var_name, var_resource, var_zone, var_units,
                          bank[1].p, scenario_label, Float64(cvar_alpha), Int(seed), Float64(scaling_factor_cost), collect(picks))
end

# Resource, zone and units of each first-stage variable, in the [x; x_line] order of capacity_variable_names.
function capacity_variable_metadata(inputs::Dict)
    names = capacity_variable_names(inputs)
    R = length(inputs["Resources"])
    zone_map, line_map = inputs["Resource zone map"], Matrix(inputs["Zone map"])
    resource = [replace(inputs["Resources"][r], r"-z\d+$" => "") for r in 1:R]
    zone = ["z$(findfirst(!=(0), zone_map[r, :]))" for r in 1:R]
    line_zone = [join(["z$z" for z in findall(!=(0), line_map[l, :])], "-") for l in 1:inputs["Number of lines"]]
    return (name = names, resource = vcat(resource, fill("Line", length(line_zone))),
            zone = vcat(zone, line_zone), units = fill("MW", length(names)))
end

function build_cut_export_input(cut_pool, MP, bank::AbstractVector, picks::AbstractVector, inputs::Dict, settings::Dict; run_id::String)
    pool_cuts, skipped = pool_optimality_cuts(cut_pool, MP, settings)
    meta = capacity_variable_metadata(inputs)
    S, F, K = length(inputs["Demand scenario probabilities"]), length(inputs["Gas price scenario probabilities"]), length(inputs["Weather scenario probabilities"])
    labels = vec(["D_$(s)_F_$(f)_W_$(k)" for s in 1:S, f in 1:F, k in 1:K])   # column-major, same as θ
    return build_cut_export_input(pool_cuts, skipped, bank, picks; run_id = run_id, var_name = meta.name, var_resource = meta.resource,
                                  var_zone = meta.zone, var_units = meta.units, scenario_label = labels,
                                  cvar_alpha = settings["Value-at-Risk percent"], seed = settings["Seed"],
                                  scaling_factor_cost = settings["Scaling factor cost"])
end

#---- Settings ------------------------------------------------------------------------------

# Hard-indexes every key so a missing setting fails loudly, like CutPool. Keeps numbers only.
function cut_export_config(settings::AbstractDict)
    return (three_parent = Int(settings["3-parent blends"]), two_parent = Int(settings["2-parent blends"]),
            sparse = Int(settings["Sparse blends"]), max_sparse_parents = Int(settings["Max sparse parents"]),
            held_out = Int(settings["Held-out blends"]), max_repair_rounds = Int(settings["Max repair rounds"]),
            clean_rounds = Int(settings["Clean rounds required"]),
            gap_tol = Float64(settings["Gap tolerance"]), near_active_tol = Float64(settings["Near-active tolerance"]),
            duplicate_tol = Float64(settings["Duplicate tolerance"]), float32_tol = Float64(settings["Float32 tolerance"]),
            validity_tol = Float64(settings["Validity tolerance"]), sparse_density = Float64(settings["Sparse density threshold"]),
            target_mb = Float64(settings["Target MB"]), hard_cap_mb = Float64(settings["Hard cap MB"]))
end

#---- 3. Prune ------------------------------------------------------------------------------

# Relative difference, safe near zero
rel_diff(a, b) = abs(a - b) / max(abs(b), 1.0)

# Random blends of `n_parents` columns of X (n × N), with Dirichlet(1,...,1) weights: normalized
# exponential draws, the same construction as random_parent_blends. Returns n × m.
function random_blends(X::AbstractMatrix, m::Int, n_parents::Int, rng::AbstractRNG)
    out = zeros(size(X, 1), m)
    for j in 1:m
        weights = randexp(rng, n_parents)
        weights ./= sum(weights)
        for (w, parent) in zip(weights, rand(rng, 1:size(X, 2), n_parents))
            out[:, j] .+= w .* view(X, :, parent)
        end
    end
    return out
end

# Sparse blends: a different number of parents, from 2 up to max_parents, for each point.
function random_sparse_blends(X::AbstractMatrix, m::Int, max_parents::Int, rng::AbstractRNG)
    m == 0 && return zeros(size(X, 1), 0)
    return reduce(hcat, [random_blends(X, 1, rand(rng, 2:max(2, max_parents)), rng) for _ in 1:m])
end

# Probe set: every solved plan plus random blends of them. density scales the blend counts.
function build_probe_blends(X::AbstractMatrix, cfg, density::Real, rng::AbstractRNG)
    count(c) = round(Int, c * density)
    return hcat(X, random_blends(X, count(cfg.three_parent), 3, rng), random_blends(X, count(cfg.two_parent), 2, rng),
                random_sparse_blends(X, count(cfg.sparse), cfg.max_sparse_parents, rng))
end

# Blends the coverage check uses: a mix of the same kinds, in the same proportions as the probes.
function build_held_out_blends(X::AbstractMatrix, cfg, rng::AbstractRNG)
    total = cfg.three_parent + cfg.two_parent + cfg.sparse
    share(c) = total == 0 ? 0 : round(Int, cfg.held_out * c / total)
    return hcat(random_blends(X, share(cfg.three_parent), 3, rng), random_blends(X, share(cfg.two_parent), 2, rng),
                random_sparse_blends(X, share(cfg.sparse), cfg.max_sparse_parents, rng))
end

# Largest cut value at each column of X, and the cut attaining it. Works in chunks of probes so
# the temporary K × chunk matrix stays small however many probes there are.
function best_cut_at(alpha::AbstractVector, beta::AbstractMatrix, X::AbstractMatrix)
    values, cuts = zeros(size(X, 2)), zeros(Int, size(X, 2))
    for cols in Iterators.partition(axes(X, 2), CUT_EXPORT_CHUNK_PROBES)
        cut_values_here = alpha .+ beta' * X[:, cols]               # K × chunk
        best = vec(argmax(cut_values_here; dims = 1))
        values[cols] = cut_values_here[best]
        cuts[cols] = getindex.(best, 1)
    end
    return values, cuts
end

# Which cuts are within near_tol (relative) of the largest cut at some column of X. near_tol = 0
# keeps only the maximizers, one per column.
function active_cut_mask(alpha::AbstractVector, beta::AbstractMatrix, X::AbstractMatrix, near_tol::Real)
    active = falses(length(alpha))
    for cols in Iterators.partition(axes(X, 2), CUT_EXPORT_CHUNK_PROBES)
        cut_values_here = alpha .+ beta' * X[:, cols]               # K × chunk
        best = vec(argmax(cut_values_here; dims = 1))
        active[getindex.(best, 1)] .= true
        if near_tol > 0
            top = cut_values_here[best]
            threshold = top .- near_tol .* max.(abs.(top), 1.0)
            active .|= vec(any(cut_values_here .>= threshold'; dims = 2))
        end
    end
    return active
end

# Indices of the cuts that survive deduplication, sorted.
#   - Identical slopes: only the largest intercept can ever be the max.
#   - Near-duplicates: normalize (α, β·xscale) to unit length, so that both parts are in $ at a
#     typical capacity, round to `tol`, and keep the largest intercept of each group.
# Dropping a near-duplicate can lower the max a little; the coverage check measures that because
# it always compares against the full archive.
function dedupe_cuts(alpha::AbstractVector, beta::AbstractMatrix, xscale::AbstractVector, tol::Real)
    best_of_group = Dict{Any, Int}()
    function keep_best!(key, k)
        current = get(best_of_group, key, 0)
        (current == 0 || alpha[k] > alpha[current]) && (best_of_group[key] = k)
    end
    for k in eachindex(alpha)
        keep_best!(beta[:, k], k)
    end
    survivors = collect(values(best_of_group))
    empty!(best_of_group)
    for k in survivors
        c = vcat(alpha[k], beta[:, k] .* xscale)
        keep_best!(round.(Int, c ./ max(norm(c), eps()) ./ tol), k)
    end
    return sort!(collect(values(best_of_group)))
end

# Prune one scenario; returns the kept cut indices (sorted) and the repair rounds used.
# Repair stops after `clean_rounds` held-out sets in a row expose no gap. One clean set isn't enough
# evidence: gaps can sit in small regions of the blend space that a single sample easily misses.
function prune_scenario(alpha::AbstractVector, beta::AbstractMatrix, probes::AbstractMatrix, held_out_sets::AbstractVector,
                        xscale::AbstractVector, near_tol::Real, gap_tol::Real, duplicate_tol::Real, clean_rounds::Int)
    kept = falses(length(alpha))

    # Near-active cuts at the probes, chosen among the deduplicated cuts
    unique_cuts = dedupe_cuts(alpha, beta, xscale, duplicate_tol)
    kept[unique_cuts[active_cut_mask(alpha[unique_cuts], beta[:, unique_cuts], probes, near_tol)]] .= true

    # Coverage check and repair: add the true maximizer wherever a fresh blend exposes a gap
    rounds, clean_in_a_row = 0, 0
    for held_out in held_out_sets
        rounds += 1
        full_values, full_cuts = best_cut_at(alpha, beta, held_out)
        pruned_values, _ = best_cut_at(alpha[kept], beta[:, kept], held_out)
        exposed = rel_diff.(pruned_values, full_values) .> gap_tol
        kept[full_cuts[exposed]] .= true
        clean_in_a_row = any(exposed) ? 0 : clean_in_a_row + 1
        clean_in_a_row >= clean_rounds && break
    end
    return findall(kept), rounds
end

# Prune every scenario, streaming through them so only one scenario's temporaries exist at a time.
# Reports the final gaps on `test_blends`, which no repair round has seen.
function prune_cut_table(table::CutTable, probes::AbstractMatrix, held_out_sets::AbstractVector, test_blends::AbstractMatrix,
                         xscale::AbstractVector, near_tol::Real, gap_tol::Real, duplicate_tol::Real, clean_rounds::Int)
    S = length(table.alpha)
    kept = Vector{Vector{Int}}(undef, S)
    rounds = zeros(Int, S)
    gaps = Float64[]
    for s in 1:S
        alpha, beta = table.alpha[s], table.beta[s]
        kept[s], rounds[s] = prune_scenario(alpha, beta, probes, held_out_sets, xscale, near_tol, gap_tol, duplicate_tol, clean_rounds)
        full_values, _ = best_cut_at(alpha, beta, test_blends)
        pruned_values, _ = best_cut_at(alpha[kept[s]], beta[:, kept[s]], test_blends)
        append!(gaps, rel_diff.(pruned_values, full_values))
    end
    return kept, (max_gap = maximum(gaps; init = 0.0), mean_gap = isempty(gaps) ? 0.0 : mean(gaps),
                  q99_gap = isempty(gaps) ? 0.0 : quantile(gaps, 0.99), share_over_tol = isempty(gaps) ? 0.0 : mean(gaps .> gap_tol),
                  max_rounds_used = maximum(rounds), n_points = length(gaps))
end

#---- 4. Write ------------------------------------------------------------------------------

# Largest relative change in LB_s at any plan if β and/or x are stored as Float32. The dashboard
# computes in Float64, so only the rounding of the stored numbers counts: values are rounded to
# Float32 and converted back before multiplying. α is always Float64.
function float32_lb_change(table::CutTable, kept, X::AbstractMatrix; round_beta::Bool, round_x::Bool)
    worst = 0.0
    rounded(A, on) = on ? Float64.(Float32.(A)) : A
    X_stored = rounded(X, round_x)
    for s in eachindex(kept)
        alpha, beta = table.alpha[s][kept[s]], table.beta[s][:, kept[s]]
        exact, _ = best_cut_at(alpha, beta, X)
        stored, _ = best_cut_at(alpha, rounded(beta, round_beta), X_stored)
        worst = max(worst, maximum(rel_diff.(stored, exact); init = 0.0))
    end
    return worst
end

# Smallest storage that keeps LB within tol at every plan, trying x before β in Float64 because x is the smaller array.
function choose_float_types(table::CutTable, kept, X::AbstractMatrix, tol::Real)
    changes = Dict{String, Float64}()
    for (round_beta, round_x) in ((true, true), (true, false), (false, true))
        change = float32_lb_change(table, kept, X; round_beta = round_beta, round_x = round_x)
        changes["beta $(round_beta ? 32 : 64), x $(round_x ? 32 : 64)"] = change
        change <= tol && return (beta = round_beta ? Float32 : Float64, x = round_x ? Float32 : Float64, changes = changes)
    end
    return (beta = Float64, x = Float64, changes = changes)
end

# Numeric array with gzip and the shuffle filter. Chunks are whole columns (cuts or plans).
function write_compressed(group, name::AbstractString, data::AbstractArray)
    if isempty(data)
        group[name] = data
        return
    end
    chunk = ndims(data) == 1 ? (min(length(data), 16384),) : (size(data, 1), min(size(data, 2), 4096))
    HDF5.write_dataset(group, name, data; chunk = chunk, shuffle = true, deflate = 6)
end

# Writes the kept cuts of every scenario. Everything the budget loop can vary is an argument.
function write_cut_h5(path::AbstractString, input::CutExportInput, kept, float_types, cfg, step, coverage)
    table = input.cuts
    n = size(input.plan_x, 1)
    K = sum(length, kept)
    scenario = Vector{Int16}(undef, K)
    alpha = Vector{Float64}(undef, K)
    iteration, source, batch = Vector{Int16}(undef, K), Vector{Int16}(undef, K), Vector{Int16}(undef, K)
    beta = Matrix{float_types.beta}(undef, n, K)
    position = 0
    for s in eachindex(kept), k in kept[s]
        position += 1
        scenario[position], alpha[position] = s, table.alpha[s][k]
        iteration[position], source[position], batch[position] = table.iteration[s][k], table.source[s][k], table.batch[s][k]
        beta[:, position] = table.beta[s][:, k]
    end
    density = K == 0 ? 0.0 : count(!iszero, beta) / length(beta)
    use_csr = density < cfg.sparse_density

    HDF5.h5open(path, "w") do f
        cuts = HDF5.create_group(f, "cuts")
        write_compressed(cuts, "scenario", scenario)
        write_compressed(cuts, "alpha", alpha)
        write_compressed(cuts, "iteration", iteration)
        write_compressed(cuts, "source", source)
        write_compressed(cuts, "batch", batch)
        if use_csr
            # Each Julia column is one cut, so the CSC arrays of the n×K matrix are the CSR arrays of K×n
            csr = HDF5.create_group(cuts, "beta_csr")
            sparse_beta = sparse(beta)
            write_compressed(csr, "indptr", Int64.(sparse_beta.colptr .- 1))
            write_compressed(csr, "indices", Int32.(sparse_beta.rowval .- 1))
            write_compressed(csr, "data", sparse_beta.nzval)
        else
            write_compressed(cuts, "beta", beta)
        end

        vars = HDF5.create_group(f, "vars")
        vars["name"], vars["resource"], vars["zone"], vars["units"] = input.var_name, input.var_resource, input.var_zone, input.var_units

        plans = HDF5.create_group(f, "plans")
        plans["id"] = input.plan_ids
        write_compressed(plans, "x", Matrix{float_types.x}(input.plan_x))

        meta = HDF5.create_group(f, "meta")
        meta["run_id"] = input.run_id
        meta["n_scenarios"] = length(input.p)
        meta["scenario_probability"] = input.p
        meta["scenario_label"] = input.scenario_label
        meta["sources"] = input.sources
        meta["cvar_alpha"] = input.cvar_alpha
        meta["alpha_units"] = "\$"
        meta["beta_units"] = "\$/MW"
        meta["theta_definition"] = "annual undiscounted operating cost: variable generation + non-served energy + capacity-reserve slack; no investment"
        meta["master_scaling_factor_cost"] = input.scaling_factor_cost   # already removed from every number in this file
        meta["feasibility_cuts"] = "none: subproblems have complete recourse"
        meta["skipped_pool_rows"] = input.skipped_rows
        meta["failed_plans"] = input.plan_ids[.!input.plan_ok]
        meta["raw_cuts_per_scenario"] = length.(table.alpha)
        meta["pruned_cuts_per_scenario"] = length.(kept)
        meta["raw_cuts_total"] = n_table_cuts(table)
        meta["pruned_cuts_total"] = K
        meta["beta_format"] = use_csr ? "csr" : "dense"
        meta["beta_density"] = density
        meta["beta_type"], meta["x_type"] = string(float_types.beta), string(float_types.x)

        pruning = HDF5.create_group(meta, "pruning")
        pruning["near_active_tol"] = step.near_active_tol
        pruning["gap_tol"] = step.gap_tol
        pruning["probe_density"] = step.density
        pruning["duplicate_tol"] = cfg.duplicate_tol
        pruning["three_parent_blends"], pruning["two_parent_blends"], pruning["sparse_blends"] = cfg.three_parent, cfg.two_parent, cfg.sparse
        pruning["max_sparse_parents"], pruning["held_out_blends"], pruning["max_repair_rounds"] = cfg.max_sparse_parents, cfg.held_out, cfg.max_repair_rounds
        pruning["clean_rounds_required"] = cfg.clean_rounds

        gap = HDF5.create_group(meta, "coverage_gap")
        gap["max"], gap["mean"], gap["q99"] = coverage.max_gap, coverage.mean_gap, coverage.q99_gap
        gap["share_over_tol"], gap["max_repair_rounds_used"], gap["n_points"] = coverage.share_over_tol, coverage.max_rounds_used, coverage.n_points
    end
    return path
end

# The cuts in a written file, as a CutTable, plus the plans. Used to validate what was really written.
function read_cut_h5(path::AbstractString)
    HDF5.h5open(path, "r") do f
        scenario = read(f["cuts/scenario"])
        alpha, iteration, source, batch = read(f["cuts/alpha"]), read(f["cuts/iteration"]), read(f["cuts/source"]), read(f["cuts/batch"])
        S = Int(read(f["meta/n_scenarios"]))
        n = length(read(f["vars/name"]))
        K = length(alpha)
        beta = if read(f["meta/beta_format"]) == "csr"
            indptr, indices, data = read(f["cuts/beta_csr/indptr"]), read(f["cuts/beta_csr/indices"]), read(f["cuts/beta_csr/data"])
            Matrix(SparseMatrixCSC(n, K, Int.(indptr) .+ 1, Int.(indices) .+ 1, Float64.(data)))
        else
            Float64.(read(f["cuts/beta"]))
        end
        rows = [findall(==(s), scenario) for s in 1:S]
        table = CutTable([alpha[r] for r in rows], [beta[:, r] for r in rows], [Int.(source[r]) for r in rows],
                         [Int.(iteration[r]) for r in rows], [Int.(batch[r]) for r in rows])
        return (cuts = table, plan_ids = read(f["plans/id"]), plan_x = Float64.(read(f["plans/x"])))
    end
end

# Stored (compressed) bytes of every dataset, by path.
function dataset_sizes(path::AbstractString)
    sizes = Pair{String, Int}[]
    function walk!(group, prefix)
        for key in keys(group)
            object = group[key]
            if object isa HDF5.Group
                walk!(object, prefix * key * "/")
            else
                push!(sizes, prefix * key => Int(HDF5.API.h5d_get_storage_size(object)))
            end
        end
    end
    HDF5.h5open(f -> walk!(f, ""), path, "r")
    return sort!(sizes; by = last, rev = true)
end

#---- 5. Validate ---------------------------------------------------------------------------

# The cuts that existed when GridFill batch b was aimed: everything from before b, in table order.
function cut_bank_before_batch(table::CutTable, b::Integer)
    masks = [table.batch[s] .< b for s in eachindex(table.alpha)]
    return CutBank([table.alpha[s][masks[s]] for s in eachindex(masks)], [table.beta[s][:, masks[s]] for s in eachindex(masks)])
end

# Recomputes GF_LB_EV and GF_LB_CVaR by repeating exactly what the sampler did for each batch:
# the same plans and cuts, the same prune_cuts with the same seeded blends, and the same clamp
# inside blend_bracket. These logged bounds use a pruned subset of that batch's cuts, so this is
# the only computation that can match them exactly.
function replay_logged_bounds(input::CutExportInput; prune_blends::Int = GRIDFILL_TUNING.prune_blends)
    picks = input.picks
    lb_ev, lb_cvar = zeros(length(picks)), zeros(length(picks))
    for b in unique(pk.batch for pk in picks)
        solved = findall(input.plan_ok .& (input.plan_batch .< b))
        X = input.plan_x[:, solved]
        index_of = Dict(input.plan_ids[j] => i for (i, j) in enumerate(solved))
        Xref = hcat(X, random_parent_blends(permutedims(X), prune_blends, Xoshiro(hash((input.seed, b)))))
        cuts = prune_cuts(cut_bank_before_batch(input.cuts, b), Xref)
        members = findall(pk -> pk.batch == b, picks)
        parents = reduce(vcat, [[index_of[name] for name in picks[i].parent_names]' for i in members])
        W = reduce(vcat, [picks[i].weights' for i in members])
        Xb = reduce(hcat, [picks[i].x for i in members])
        br = blend_bracket(cuts, input.plan_θ[:, solved], zeros(length(solved)), input.p, input.cvar_alpha, parents, W, Xb)
        lb_ev[members], lb_cvar[members] = br.lb_ev, br.lb_cvar
    end
    return (lb_ev = lb_ev, lb_cvar = lb_cvar)
end

# EV and CVaR floors at each GridFill pick, from the cuts of `table` older than the pick's batch.
function pick_floors(table::CutTable, input::CutExportInput)
    picks = input.picks
    pn = input.p ./ sum(input.p)
    ev, cvar = zeros(length(picks)), zeros(length(picks))
    for b in unique(pk.batch for pk in picks)
        members = findall(pk -> pk.batch == b, picks)
        L = cut_values(cut_bank_before_batch(table, b), reduce(hcat, [picks[i].x for i in members]))
        ev[members] = vec(pn' * L)
        cvar[members] = cvar_columns(L, pn, input.cvar_alpha)
    end
    return (ev = ev, cvar = cvar)
end

max_mean_rel(a, b) = isempty(a) ? (0.0, 0.0) : (maximum(rel_diff.(a, b)), mean(rel_diff.(a, b)))

# Step 3 of the spec. Reads the written file back, so it checks what the dashboard will load.
function validate_export(path::AbstractString, input::CutExportInput, cfg; prune_blends::Int = GRIDFILL_TUNING.prune_blends)
    written = read_cut_h5(path)
    lines = String[]
    say(text) = push!(lines, text)
    say("Cut archive report: $(input.run_id)")
    say(@sprintf("Cuts: %d raw -> %d pruned, %d scenarios, %d variables, %d plans", n_table_cuts(input.cuts), n_table_cuts(written.cuts), length(input.p), size(input.plan_x, 1), length(input.plan_ids)))
    say("Feasibility cuts: none (complete recourse). Pool rows with no alpha term, skipped: $(input.skipped_rows).")
    !all(input.plan_ok) && say("Plans with failed subproblems, no checks and no cuts: $(join(input.plan_ids[.!input.plan_ok], ", "))")

    # 1. Logged GridFill bounds
    replay_diff = (ev = (0.0, 0.0), cvar = (0.0, 0.0))
    if isempty(input.picks)
        say("1. Logged bounds: GridFill didn't run, nothing to reproduce.")
    else
        logged_ev, logged_cvar = [pk.lb_ev for pk in input.picks], [pk.lb_cvar for pk in input.picks]
        replay = replay_logged_bounds(input; prune_blends = prune_blends)
        replay_diff = (ev = max_mean_rel(replay.lb_ev, logged_ev), cvar = max_mean_rel(replay.lb_cvar, logged_cvar))
        full, pruned = pick_floors(input.cuts, input), pick_floors(written.cuts, input)
        say(@sprintf("1. Logged bounds at %d GridFill picks, relative difference (max, mean):", length(input.picks)))
        say(@sprintf("   replay of the sampler vs logged:  EV %.2e, %.2e   CVaR %.2e, %.2e   (should be ~0)", replay_diff.ev..., replay_diff.cvar...))
        say(@sprintf("   full archive at that batch vs logged:  EV %.2e, %.2e   CVaR %.2e, %.2e   (>= logged: the sampler used a pruned subset and a clamp)", max_mean_rel(full.ev, logged_ev)..., max_mean_rel(full.cvar, logged_cvar)...))
        say(@sprintf("   pruned file at that batch vs full archive:  EV %.2e, %.2e   CVaR %.2e, %.2e   (overstates the loss: pruning targets the final archive, not this earlier state)", max_mean_rel(pruned.ev, full.ev)..., max_mean_rel(pruned.cvar, full.cvar)...))
    end

    # 2 and 3. Validity and tightness at the plans, from the plans' own subproblem solves
    ok = findall(input.plan_ok)
    X, θ = input.plan_x[:, ok], input.plan_θ[:, ok]
    report_excess(label, L) = begin
        excess = (L .- θ) ./ max.(abs.(θ), 1.0)
        n_over = count(>(cfg.validity_tol), excess)
        i = argmax(excess)
        say(@sprintf("2. Validity, %s: %d of %d (plan, scenario) pairs exceed Op_s by more than %.1e; worst %.2e relative, \$%.3g (%s, %s)",
                     label, n_over, length(excess), cfg.validity_tol, excess[i], (L .- θ)[i], input.plan_ids[ok[i[2]]], input.scenario_label[i[1]]))
        return n_over
    end
    violations_raw = report_excess("full archive", cut_values(CutBank(input.cuts.alpha, input.cuts.beta), X))
    L_file = cut_values(CutBank(written.cuts.alpha, written.cuts.beta), written.plan_x[:, ok])
    violations_file = report_excess("written file", L_file)
    slack = (θ .- L_file) ./ max.(abs.(θ), 1.0)
    say(@sprintf("3. Tightness of the written file, (Op_s - LB_s) / Op_s: max %.2e, mean %.2e; share below 1e-6: %.3f", maximum(abs, slack), mean(abs, slack), mean(abs.(slack) .< 1e-6)))

    # 4. Sizes
    sizes = dataset_sizes(path)
    say(@sprintf("4. File size: %.2f MB. Stored size per dataset:", filesize(path) / 1e6))
    for (dataset, bytes) in sizes
        bytes >= 1024 && say(@sprintf("   %-28s %9.1f KB", dataset, bytes / 1024))
    end
    return (lines = lines, replay_diff = replay_diff, violations_raw = violations_raw, violations_file = violations_file,
            max_slack = maximum(abs, slack), written = written)
end

#---- The whole export ----------------------------------------------------------------------

# Shrink steps, used in order until the file fits the target. A larger near-active tolerance
# keeps more cuts, so shrinking lowers it. The probes come next, and the gap tolerance last
# because it is the only step that accepts a worse bound.
function budget_steps(cfg)
    steps = [(near_active_tol = cfg.near_active_tol, density = 1.0, gap_tol = cfg.gap_tol),
             (near_active_tol = 0.0, density = 1.0, gap_tol = cfg.gap_tol),
             (near_active_tol = 0.0, density = 0.5, gap_tol = cfg.gap_tol),
             (near_active_tol = 0.0, density = 0.25, gap_tol = cfg.gap_tol),
             (near_active_tol = 0.0, density = 0.25, gap_tol = 10 * cfg.gap_tol)]
    return unique(steps)
end

format_budget_row(step, K, mb, coverage) = @sprintf("near-active %.1e, probe density %.2f, gap tol %.1e: %d cuts, %.2f MB, gap at %d unseen points: max %.2e, mean %.2e, share over tol %.5f%%",
                                                    step.near_active_tol, step.density, step.gap_tol, K, mb, coverage.n_points, coverage.max_gap, coverage.mean_gap, 100 * coverage.share_over_tol)

# Tries budget_steps until the file fits the target. Returns the rows of the size-vs-gap tradeoff;
# `output_path` holds the chosen file. With nothing under the hard cap, it holds no file and this errors.
function write_within_budget(output_path::AbstractString, input::CutExportInput, cfg)
    usable = findall(input.plan_ok)
    isempty(usable) && error("Cut export: every plan has failed subproblems.")
    X = input.plan_x[:, usable]
    xscale = [max(v, 1.0) for v in vec(maximum(abs, X; dims = 2))]   # typical capacity of each variable, MW

    # Held-out and final test blends stay the same across steps, so the steps can be compared
    held_rng = Xoshiro(input.seed + 1)
    held_out_sets = [build_held_out_blends(X, cfg, held_rng) for _ in 1:cfg.max_repair_rounds]
    test_blends = build_held_out_blends(X, cfg, held_rng)

    candidate = output_path * ".tmp"
    rows = String[]
    best_mb = Inf
    for step in budget_steps(cfg)
        probes = build_probe_blends(X, cfg, step.density, Xoshiro(input.seed))
        kept, coverage = prune_cut_table(input.cuts, probes, held_out_sets, test_blends, xscale, step.near_active_tol, step.gap_tol, cfg.duplicate_tol, cfg.clean_rounds)
        float_types = choose_float_types(input.cuts, kept, X, cfg.float32_tol)
        for (combo, change) in sort!(collect(float_types.changes))
            change > cfg.float32_tol && @warn(@sprintf("Cut export: storing %s changes LB by up to %.2e at the plans, over the %.1e tolerance; using more precision.", combo, change, cfg.float32_tol))
        end
        write_cut_h5(candidate, input, kept, float_types, cfg, step, coverage)
        mb = filesize(candidate) / 1e6
        push!(rows, format_budget_row(step, sum(length, kept), mb, coverage))
        @info("Cut export budget step: " * rows[end])
        if mb < best_mb
            best_mb = mb
            mv(candidate, output_path; force = true)
        end
        mb <= cfg.target_mb && break
    end
    rm(candidate; force = true)
    if best_mb > cfg.hard_cap_mb
        rm(output_path; force = true)
        error("Cut export: smallest file is $(round(best_mb; digits = 2)) MB, over the $(cfg.hard_cap_mb) MB hard cap. No file written. Tradeoff:\n" * join(rows, "\n"))
    end
    best_mb > cfg.target_mb && @warn("Cut export: smallest file is $(round(best_mb; digits = 2)) MB, over the $(cfg.target_mb) MB target but under the hard cap.")
    return rows
end

# Runs the whole export and writes <file_stem>.h5 in `folder`. The raw copy is removed on success
# and left in place on failure; rerun with export_cut_archive_from_raw.
function export_cut_archive(input::CutExportInput, folder::AbstractString, settings::AbstractDict;
                            file_stem::AbstractString = "cut_archive", prune_blends::Int = GRIDFILL_TUNING.prune_blends)
    cfg = cut_export_config(settings)
    mkpath(folder)
    raw_path = joinpath(folder, file_stem * "_raw.jls")
    output_path = joinpath(folder, file_stem * ".h5")
    try
        serialize(raw_path, input)
    catch e
        @warn("Cut export: couldn't save the raw copy; continuing without one", exception = (e, catch_backtrace()))
    end
    prune_seconds = @elapsed tradeoff = write_within_budget(output_path, input, cfg)
    validate_seconds = @elapsed report = validate_export(output_path, input, cfg; prune_blends = prune_blends)
    @info(@sprintf("Cut export: pruning and writing %.1f s, validation %.1f s", prune_seconds, validate_seconds))
    HDF5.h5open(output_path, "r+") do f
        f["meta/budget_tradeoff"] = tradeoff
        f["meta/validation_report"] = report.lines
    end
    println(join(report.lines, "\n"))
    println("Size vs gap at each budget step:\n" * join(tradeoff, "\n"))
    rm(raw_path; force = true)
    return (path = output_path, tradeoff = tradeoff, report = report)
end

# Prune again from the raw copy a failed export left behind, without the model.
export_cut_archive_from_raw(raw_path::AbstractString, folder::AbstractString, settings::AbstractDict; kwargs...) =
    export_cut_archive(deserialize(raw_path), folder, settings; kwargs...)
