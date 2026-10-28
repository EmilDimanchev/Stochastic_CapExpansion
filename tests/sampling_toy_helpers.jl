using LinearAlgebra, Random

"""
Synthetic data shared by the solver-free sampler tests (boundary_sampling_tests.jl,
gridfill_tests.jl): each scenario's op cost is a convex piecewise-linear function of capacity,
θ_s(x) = max_l (a_sl + c_slᵀ x), so its exact gradient (the active piece's c) is known and
the bounds can be checked against the true value of any blend.
"""

#---- Synthetic convex op-cost model -------------------------------------------------------
struct ToyCost
    a::Array{Float64, 2}   # S × L intercepts
    c::Array{Float64, 3}   # S × L × d slopes
end
function toy_cost(rng, S, L, d)
    return ToyCost(100 .+ 50 .* rand(rng, S, L), -rand(rng, S, L, d))
end
function toy_eval(tc::ToyCost, x)
    S, L, d = size(tc.c)
    θ, β = zeros(S), zeros(d, S)
    for s in 1:S
        vals = [tc.a[s, l] + dot(tc.c[s, l, :], x) for l in 1:L]
        l = argmax(vals)
        θ[s], β[:, s] = vals[l], tc.c[s, l, :]
    end
    return θ, β
end

# A bank with every solve family. x is random in [0, 20]^d; θ and β come from the toy model.
function toy_bank(rng; d = 5, S = 20, n_vectors = 8, n_budgets = 3)
    tc = toy_cost(rng, S, 4, d)
    p = fill(1 / S, S)
    names = ["Risk_Weight_0.0", "Risk_Weight_1.0", "Risk_Weight_0.5"]
    for k in 1:n_vectors, b in 1:n_budgets
        push!(names, "Random_$(k)_Budget_$(b)")
        if b == 1
            append!(names, ["Random_$(k)_Budget_1_System_Expected_Cap_$j" for j in 1:2])
            append!(names, ["Random_$(k)_Budget_1_CVaR_Cap_$j" for j in 1:2])
        end
    end
    bank = FinalSolve[]
    for n in names
        x = 20 .* rand(rng, d)
        θ, β = toy_eval(tc, x)
        push!(bank, FinalSolve(n, "", x, θ, β, p, sum(x), true))
    end
    assign_regions!(bank)
    return bank, tc
end
