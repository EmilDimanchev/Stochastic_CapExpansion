using LinearAlgebra, Random

"""
Synthetic data for the solver-free sampler tests (gridfill_tests.jl): each scenario's op cost
is a convex piecewise-linear function of capacity, θ_s(x) = max_l (a_sl + c_slᵀ x), so its exact
gradient (the active piece's c) is known and the bounds can be checked against the true value
of any blend.
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
