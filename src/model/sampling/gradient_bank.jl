using LinearAlgebra, Random, Statistics, Dates, DataFrames, CSV, YAML

#==========================================================================================

Gradient bank: one record per solved capacity portfolio

Every final MGA solve (and every evaluated GridFill sample) is stored as a FinalSolve. It
holds everything GridFill needs about that portfolio:

  x  - capacity vector [resource capacities; line capacities], MW, length d. This is the
       order of the master's [x; x_line] variables and of the Summary capacity columns.
  θ  - per-scenario operating cost at x, $, length S. Scenarios are flattened (s,f,k)
       column-major, i.e. the same order compute_cvar uses.
  β  - per-scenario capacity gradient ∂θ_s/∂x, $/MW, size d × S. These are the Benders
       cut coefficients: the duals of the subproblem capacity Parameters
       (write_subproblem_results: "capacity dual" and "line dual"). They are usually <= 0.
       An empty d × 0 matrix means "no gradients available".
  p  - scenario probabilities (P_s[s]*P_f[f]*P_k[k]), length S.
  θb, βb - the same per subproblem (scenario × time block, see subproblem_blocks): length S·B
       and d × S·B, flattened column-major over (s,f,k,b), so entry ω + (b−1)·S is block b of
       scenario ω. θ and β are their sums over blocks. The cuts GridFill and the cut export use
       are these per-block tangents; with one block per scenario θb = θ and βb = β.

Important: θ and β are always taken from subproblems solved AT x. The Summary columns
EV_OpCost/CVaR_OpCost come from the last Benders iterate's estimate, which can belong to a
slightly different capacity (the converged MP is re-solved without regularization before
reporting), so they are never used for GridFill's brackets.

==========================================================================================#

mutable struct FinalSolve
    name::String
    x::Vector{Float64}
    θ::Vector{Float64}
    β::Matrix{Float64}
    p::Vector{Float64}
    inv_cost::Float64
    ok::Bool                  # false if any scenario failed to solve to optimality
    θb::Vector{Float64}       # per subproblem op cost, length S·B
    βb::Matrix{Float64}       # per subproblem gradient, d × S·B
end

# One block per scenario: the subproblem terms are the scenario terms
FinalSolve(name, x, θ, β, p, inv_cost, ok) = FinalSolve(name, x, θ, β, p, inv_cost, ok, θ, β)

has_gradients(fs::FinalSolve) = size(fs.βb, 2) == length(fs.θb) && !isempty(fs.θb)

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
    ok || @warn("Solve $name has subproblems that did not solve to optimality; it will not be used as a parent.")
    # Per subproblem, scenario index fastest (see the header)
    blocks = [(ω, b) for ω in eachindex(outs), b in eachindex(outs[1]["Block SP objective"])]
    θb = vec([outs[ω]["Block SP objective"][b] for (ω, b) in blocks])
    βb = reduce(hcat, vec([vcat(outs[ω]["Block capacity dual"][:, b], outs[ω]["Block line dual"][:, b]) for (ω, b) in blocks]))
    return FinalSolve(name, copy(x), θ, Matrix{Float64}(β), p, inv_cost, ok, θb, Matrix{Float64}(βb))
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

# Optional persistence (Sampling setting "Write gradients"): one row per scenario, with
# θ_s and the d gradient entries.
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
