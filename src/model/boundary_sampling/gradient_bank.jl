using LinearAlgebra, Random, Statistics, Dates, DataFrames, CSV, YAML

#==========================================================================================

Gradient bank: one record per solved capacity portfolio

Every final MGA solve (and, later, every evaluated interior sample) is stored as a
FinalSolve. It holds everything the boundary sampler needs about that portfolio:

  x  - capacity vector [resource capacities; line capacities], MW, length d. This is the
       order of the master's [x; x_line] variables and of the Summary capacity columns.
  θ  - per-scenario operating cost at x, $, length S. Scenarios are flattened (s,f,k)
       column-major, i.e. the same order compute_cvar uses.
  β  - per-scenario capacity gradient ∂θ_s/∂x, $/MW, size d × S. These are the Benders
       cut coefficients: the duals of the subproblem capacity Parameters
       (write_subproblem_results: "capacity dual" and "line dual"). They are usually <= 0.
       An empty d × 0 matrix means "no gradients available" (lengthscale fallback).
  p  - scenario probabilities (P_s[s]*P_f[f]*P_k[k]), length S.

Important: θ and β are always taken from subproblems solved AT x. The Summary columns
EV_OpCost/CVaR_OpCost come from the last Benders iterate's estimate, which can belong to a
slightly different capacity (the converged MP is re-solved without regularization before
reporting), so they are never used for the sampler's bounds.

==========================================================================================#

mutable struct FinalSolve
    name::String
    region::String            # assigned by assign_regions! (or "interior_sample")
    x::Vector{Float64}
    θ::Vector{Float64}
    β::Matrix{Float64}
    p::Vector{Float64}
    inv_cost::Float64
    ok::Bool                  # false if any scenario failed to solve to optimality
end

has_gradients(fs::FinalSolve) = size(fs.β, 2) == length(fs.θ) && !isempty(fs.θ)

# Expected and CVaR operating cost of a solve, from its own per-scenario costs
solve_ev(fs::FinalSolve) = dot(fs.p, fs.θ) / sum(fs.p)
solve_cvar(fs::FinalSolve, α::Real) = dot(cvar_tail_weights(fs.θ, fs.p, α), fs.θ) / α

# Flattened scenario probabilities, same (s,f,k) order as the flattened θ
function scenario_probabilities(inputs::Dict)
    P_s = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]
    return vec([P_s[s]*P_f[f]*P_k[k] for s in eachindex(P_s), f in eachindex(P_f), k in eachindex(P_k)])
end

# Names of the d capacity variables, matching the Summary capacity columns
capacity_variable_names(inputs::Dict) = vcat(inputs["Resources"], ["Line " * string(l) for l in 1:inputs["Number of lines"]])

# Investment cost is linear in capacity: the same expression as the master's inv_cost and
# inv_cost_by_zone (master_planning.jl), so no master solve is needed for a fixed portfolio.
function investment_cost(x::AbstractVector, inputs::Dict)
    R = length(inputs["Resources"])
    cap, line = x[1:R], x[R+1:end]
    total = sum(cap .* inputs["Investment costs"]) + sum(line .* inputs["CAPEX per MW"]; init = 0.0)
    zone_map = inputs["Resource zone map"]
    by_zone = [sum(cap[r]*inputs["Investment costs"][r]*zone_map[r,z] for r in 1:R) for z in 1:inputs["Number of zones"]]
    return total, by_zone
end

# Build a FinalSolve from the output array of run_all_subproblems at capacity x
function final_solve_from_sp_outputs(name::String, x::Vector{Float64}, sp_outputs, p::Vector{Float64}, inv_cost::Float64)
    outs = vec(sp_outputs)   # (s,f,k) column-major, same as p
    θ = [o["SP objective"] for o in outs]
    β = reduce(hcat, [vcat(o["capacity dual"], o["line dual"]) for o in outs])
    ok = all(o["coeff"] == 1 for o in outs)
    ok || @warn("Solve $name has subproblems that did not solve to optimality; it will not be used in the metric or as a parent.")
    return FinalSolve(name, "", copy(x), θ, Matrix{Float64}(β), p, inv_cost, ok)
end

# Called on a Benders output dict just before release_heavy_payload! deletes "SPs".
# output["SPs"] holds the subproblems re-solved at output["MP"]["Capacity"] at convergence.
function record_final_solve!(bank::AbstractVector, name::String, output::Dict, inputs::Dict)
    mp = output["MP"]
    x = vcat(mp["Capacity"], mp["Line expansion"])
    fs = final_solve_from_sp_outputs(name, x, output["SPs"], scenario_probabilities(inputs), Float64(mp["Inv_cost"]))
    push!(bank, fs)
    return fs
end

# Optional persistence (boundary setting "Write gradients"): one row per scenario, with
# θ_s and the d gradient entries. Lets later sessions refit the metric without re-solving.
function write_gradient_file(folder::String, fs::FinalSolve, var_names::Vector{String})
    mkpath(folder)
    df = DataFrame(Scenario = 1:length(fs.θ), Probability = fs.p, OpCost = fs.θ)
    if has_gradients(fs)
        for (i, v) in enumerate(var_names)
            df[!, "dOpCost_d_" * v] = fs.β[i, :]
        end
    end
    CSV.write(joinpath(folder, "Gradients.csv"), df)
end

#------------------------------------------------------------------------------------------
# Region labels (from the run-name conventions in run_stochastic_exploration_risk_pareto)
#   Random_<k>_Budget_<b>_System_Expected_Cap_<j>  -> ev_wing
#   Random_<k>_Budget_<b>_CVaR_Cap_<j>             -> cvar_wing
#   Risk_Weight_<w>                                 -> anchor
#   Random_<k>_Budget_<b> with b = loosest (min b)  -> outer_shell
#   any other chord solve                           -> interior
# A chord run without a Budget tag (no tightening) was solved at the loosest budget, so it
# counts as Budget 1.
#------------------------------------------------------------------------------------------

budget_level(name::AbstractString) = (m = match(r"_Budget_(\d+)", name); m === nothing ? 1 : parse(Int, m.captures[1]))

function assign_regions!(bank::AbstractVector)
    is_chord(fs) = !(occursin("System_Expected_Cap", fs.name) || occursin("CVaR_Cap", fs.name) || startswith(fs.name, "Risk_Weight") || fs.region == "interior_sample")
    chord_budgets = [budget_level(fs.name) for fs in bank if is_chord(fs)]
    loosest = isempty(chord_budgets) ? 1 : minimum(chord_budgets)
    for fs in bank
        fs.region == "interior_sample" && continue
        fs.region = occursin("System_Expected_Cap", fs.name) ? "ev_wing" :
                    occursin("CVaR_Cap", fs.name) ? "cvar_wing" :
                    startswith(fs.name, "Risk_Weight") ? "anchor" :
                    budget_level(fs.name) == loosest ? "outer_shell" : "interior"
    end
    return bank
end

#------------------------------------------------------------------------------------------
# Offline path: rebuild the bank from a finished run's Summary CSV (e.g. run 13), then
# recover the gradients by re-solving the subproblems at each final capacity.
#------------------------------------------------------------------------------------------

# Final solve of each run = its highest _iter_ row (rows without _iter_ count as iterate 0).
# θ comes from the per-scenario OperatingCost_D_s_F_f_W_k columns, which are written from the
# subproblems solved at the row's capacity.
function load_finals_from_summary(csv_path::String, inputs::Dict)
    df = CSV.read(csv_path, DataFrame)
    var_names = capacity_variable_names(inputs)
    S, F, K = length(inputs["Demand scenario probabilities"]), length(inputs["Gas price scenario probabilities"]), length(inputs["Weather scenario probabilities"])
    p = scenario_probabilities(inputs)
    op_cols = vec(["OperatingCost_D_$(s)_F_$(f)_W_$(k)" for s in 1:S, f in 1:F, k in 1:K])

    base_name(idx) = replace(string(idx), r"_iter_\d+$" => "")
    iter_of(idx) = (m = match(r"_iter_(\d+)$", string(idx)); m === nothing ? 0 : parse(Int, m.captures[1]))

    best = Dict{String, Int}()   # solve name -> row of its highest iterate
    order = String[]             # keep first-appearance order of solve names
    for (row, idx) in enumerate(df.Index)
        b = base_name(idx)
        if !haskey(best, b)
            push!(order, b)
            best[b] = row
        elseif iter_of(idx) > iter_of(df.Index[best[b]])
            best[b] = row
        end
    end

    # Investment cost is recomputed from x (linear, exact) so it matches what the interior
    # evaluation computes; the largest difference to the CSV column is logged as a sanity check.
    bank = FinalSolve[]
    max_inv_diff = 0.0
    for b in order
        row = best[b]
        x = Float64[df[row, v] for v in var_names]
        θ = Float64[df[row, c] for c in op_cols]
        inv, _ = investment_cost(x, inputs)
        max_inv_diff = max(max_inv_diff, abs(inv - df[row, "Investment_Cost"]) / abs(inv))
        push!(bank, FinalSolve(b, "", x, θ, zeros(length(x), 0), p, inv, true))
    end
    @info("Loaded $(length(bank)) final solves from $csv_path (max relative investment-cost difference to the CSV: $(max_inv_diff))")
    assign_regions!(bank)
    return bank
end

# Re-solve every scenario at each solve's capacity and store θ and β (replacing the CSV θ).
# Logs how far the re-solved θ is from the stored one (the LP optimum is unique, so the two
# should agree to solver tolerance even if a different solver produced the CSV).
function recompute_gradients!(bank::AbstractVector, SPs, inputs::Dict, settings::Dict)
    R = length(inputs["Resources"])
    max_rel = 0.0
    for (i, fs) in enumerate(bank)
        sp = run_all_subproblems(SPs, inputs, settings, fs.x[1:R], fs.x[R+1:end]; minimal_payload = true)
        new = final_solve_from_sp_outputs(fs.name, fs.x, sp, fs.p, fs.inv_cost)
        if !isempty(fs.θ)
            max_rel = max(max_rel, maximum(abs.(new.θ .- fs.θ) ./ max.(abs.(fs.θ), 1.0)))
        end
        fs.θ, fs.β, fs.ok = new.θ, new.β, new.ok
        i % 50 == 0 && @info("Recomputed gradients for $i / $(length(bank)) solves")
    end
    @info("Max relative difference between re-solved and stored per-scenario op cost: $(max_rel)")
    return max_rel
end
