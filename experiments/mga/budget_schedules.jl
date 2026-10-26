#===============================

Helpers for main pareto + exploration run



================================#


function compute_budget_floor(extreme_values::Dict, transform_cvar::Float64, transform_sys::Float64, floor_offset::Float64)
    transformed_vals = [transform_sys*ev["System Expected"] + transform_cvar*ev["CVaR"] for ev in values(extreme_values)]
    isempty(transformed_vals) && error("Need at least one risk-pareto anchor point to compute a budget floor.")
    return minimum(transformed_vals) + floor_offset
end

function compute_budget_tighten_schedule(budget_start::Float64, budget_floor::Float64, step_size::Float64)
    step_size <= 0 && error("Budget tighten step size must be positive.")
    budget_floor > budget_start && error("Budget floor ($budget_floor) must be <= starting budget ($budget_start) — check budget_offset vs. the anchor points' own transformed costs.")
    n_steps = ceil(Int, (budget_start - budget_floor) / step_size)
    return [max(budget_floor, budget_start - step_size*k) for k in 0:n_steps]
end


#-------------------------------------------------------------------------------------------
# Helper 1: piecewise-linear interpolation (Base Julia has no built-in interp function).
# Given sorted breakpoints xs with values ys, return the linearly interpolated y at x.
# Values of x outside [xs[1], xs[end]] are clamped to the end values.
#-------------------------------------------------------------------------------------------
function linear_interp(xs::Vector{Float64}, ys::Vector{Float64}, x::Real)
    # Clamp to the ends of the curve
    if x <= xs[1]
        return ys[1]
    elseif x >= xs[end]
        return ys[end]
    end

    # Index i of the last breakpoint with xs[i] ≤ x, so x lies in [xs[i], xs[i+1])
    i = searchsortedlast(xs, x)

    # Guard against a zero-length segment (two breakpoints at the same x)
    if xs[i+1] == xs[i]
        return ys[i]
    end

    # Fraction of the way from xs[i] to xs[i+1], then blend the two y values
    t = (x - xs[i]) / (xs[i+1] - xs[i])
    return ys[i] + t * (ys[i+1] - ys[i])
end


#-------------------------------------------------------------------------------------------
# Helper 2: lens end. Walk the polyline vertices in the given order (from one anchor inward)
# and return the arc position of the first point where T ≤ B. Within the crossing segment
# T is linear, so the exact crossing point is found by linear interpolation.
# Returns `nothing` if no point on the polyline satisfies T ≤ B.
#-------------------------------------------------------------------------------------------
function lens_end(B::Real, order::Vector{Int}, Tv::Vector{Float64}, arc::Vector{Float64})
    for k in 1:length(order)-1
        a = order[k]      # current vertex
        b = order[k+1]    # next vertex inward

        # The current vertex already satisfies the budget: the lens starts here
        if Tv[a] <= B
            return arc[a]
        end

        # The budget line is crossed inside segment a -> b: solve for the crossing point
        if Tv[b] <= B
            t = (Tv[a] - B) / (Tv[a] - Tv[b])
            return arc[a] + t * (arc[b] - arc[a])
        end
    end

    # Only the final vertex is left to check
    if Tv[order[end]] <= B
        return arc[order[end]]
    end
    return nothing
end


#-------------------------------------------------------------------------------------------
# Helper 3: cap values for one wing. Places n evenly spaced caps strictly between the lens
# end (tight_end) and the central edge (edge), excluding both endpoints, so every cap lies
# inside the outer wing (a cap at the central edge is fixed across budget levels and often
# does not bind once the budget tightens). n is limited by n_max and by the minimum spacing
# min_gap_frac * span. Caps are returned loose -> tight, which is the order to solve them in.
#-------------------------------------------------------------------------------------------
function wing_caps(tight_end::Real, edge::Real, span::Real, n_max::Int, min_gap_frac::Real)
    width = edge - tight_end

    # How many caps fit while keeping adjacent caps (and the endpoints) at least
    # min_gap_frac * span apart: n caps split the wing into n + 1 gaps
    n = min(n_max, Int(floor(width / (min_gap_frac * span))) - 1)
    if n <= 0
        return Float64[]
    end

    # i = n gives the loosest cap (one gap inside the central edge); i = 1 the tightest
    caps = Float64[]
    for i in n:-1:1
        push!(caps, tight_end + width * i / (n + 1))
    end
    return caps
end


#-------------------------------------------------------------------------------------------
# Main function.
#
# Inputs
# - extreme_values: Dict(risk_weight => Dict("System Expected" => …, "CVaR" => …)) for ALL
#   risk-weight solves. Intermediate weights are required to locate the lens ends and the
#   tangency point. The key convention does not matter: the EV-optimal solve is the one with
#   the lowest System Expected.
# - budgets: transformed budget levels, in the same units as k_ev*SystemExpected + k_cvar*CVaR.
# - k_ev, k_cvar: weights of the transformed budget.
# - n_max: maximum number of caps per wing per level.
# - min_gap_frac: minimum spacing between adjacent caps, as a fraction of that metric's span
#   (0.03 ≈ 2.5× the Benders overshoot in the reference run).
# - edge_pullback: moves each central edge from the tangency toward its anchor, as a fraction
#   of that metric's span. 0.0 = edges exactly at the tangency.
#
# Returns Dict(B => Dict("System_Expected" => caps, "CVaR" => caps)). Iterate over `budgets`
# (not over the Dict) to keep level order.
#
# Feasibility: every polyline point is a convex combination of solved portfolios, so by
# convexity of the operational problem some feasible portfolio is at least as good on both
# metrics. Each capped region {T ≤ B, cap} therefore contains a feasible point.
#-------------------------------------------------------------------------------------------
function compute_budget_schedule_wings(extreme_values::Dict, budgets::AbstractVector;
                                       k_ev::Real, k_cvar::Real,
                                       n_max::Int = 3, min_gap_frac::Real = 0.03,
                                       edge_pullback::Real = 0.05)

    # --- Block 1: collect the risk-weight solves as (System Expected, CVaR) pairs and sort
    #     them by System Expected, so the polyline runs EV-optimal -> CVaR-optimal.
    pts = Tuple{Float64,Float64}[]
    for v in values(extreme_values)
        push!(pts, (Float64(v["System Expected"]), Float64(v["CVaR"])))
    end
    pts = sort(pts; by = first)
    if length(pts) < 3
        error("need intermediate risk-weight solves to locate lens ends and the tangency point")
    end
    ev = [first(p) for p in pts]   # System Expected at each vertex (increasing)
    cv = [last(p) for p in pts]    # CVaR at each vertex (decreasing along a proper frontier)
    n_pts = length(pts)

    # --- Block 2: spans of each metric between the two extreme solves. These normalize the
    #     axes for arc length and set the minimum cap spacing.
    ev_span = ev[end] - ev[1]
    cv_span = cv[1] - cv[end]
    if !(ev_span > 0 && cv_span > 0)
        error("extreme solves must trade off EV against CVaR")
    end

    # --- Block 3: normalized arc length along the polyline (0 at EV-optimal, 1 at CVaR-optimal).
    #     Each segment length is measured with both axes scaled by their spans.
    seg = Float64[]
    for i in 1:n_pts-1
        push!(seg, hypot((ev[i+1] - ev[i]) / ev_span, (cv[i+1] - cv[i]) / cv_span))
    end
    arc = vcat(0.0, cumsum(seg)) ./ sum(seg)

    # --- Block 4: transformed budget value T at each polyline vertex.
    Tv = [k_ev * ev[i] + k_cvar * cv[i] for i in 1:n_pts]

    # --- Block 5: central edges at the tangency point (vertex with the lowest T), optionally
    #     pulled back toward each anchor. These are the same for every budget level.
    i_tan = argmin(Tv)
    ev_central_edge = ev[i_tan] - edge_pullback * ev_span   # EV-wing caps stop here
    cv_central_edge = cv[i_tan] - edge_pullback * cv_span   # CVaR-wing caps stop here

    # --- Block 6: for each budget level, find both lens ends and place the wing caps.
    order_from_ev   = collect(1:n_pts)       # walk inward from the EV-optimal anchor
    order_from_cvar = collect(n_pts:-1:1)    # walk inward from the CVaR-optimal anchor
    schedule = Dict{Float64,Dict{String,Vector{Float64}}}()

    for B in budgets
        s_ev = lens_end(B, order_from_ev, Tv, arc)
        s_cv = lens_end(B, order_from_cvar, Tv, arc)

        # EV wing: caps on System Expected, only if the lens end lies outside the central edge
        caps_ev = Float64[]
        if !isnothing(s_ev)
            ev_end = linear_interp(arc, ev, s_ev)
            if ev_end < ev_central_edge
                caps_ev = wing_caps(ev_end, ev_central_edge, ev_span, n_max, min_gap_frac)
            end
        end

        # CVaR wing: caps on CVaR, only if the lens end lies outside the central edge
        caps_cv = Float64[]
        if !isnothing(s_cv)
            cv_end = linear_interp(arc, cv, s_cv)
            if cv_end < cv_central_edge
                caps_cv = wing_caps(cv_end, cv_central_edge, cv_span, n_max, min_gap_frac)
            end
        end

        schedule[Float64(B)] = Dict("System_Expected" => caps_ev, "CVaR" => caps_cv)
    end
    return schedule
end
