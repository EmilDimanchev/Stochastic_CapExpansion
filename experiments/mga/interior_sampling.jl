#======================

Interior sampling functions


=======================#

function make_mgca_problem(points)
    rows, cols = size(points)
    points = max.(points, 0.0) # ensure all points are non-negative for convex combination

    model = Model(Ipopt.Optimizer)
    @variable(model, l[1:rows] >= 0)
    @constraint(model, c_sum, sum(l[i] for i in 1:rows) == 1)
    @variable(model, x[1:cols])
    @constraint(model, x == points' * l)
    @constraint(model, c_max_l, l .<= 0.95) # do not replicate individual points
    @objective(model, Min, 0)
    set_silent(model)
    return model
end

function sample_interior(points, num_samples, settings)
    return sample_interior_distributed(points, num_samples, settings)
end

# takes in points, cost_df, num samples and settings
# transforms cost_df into a 2d array of inv+ev and inv+cvar columns, with each row being a point in order of solve
# sends those to sample_interior_distributed which will return the evaluated mps

function sample_interior_delauney(points, cost_df::DataFrame, num_samples, settings, budget::Float64, transform_sys::Float64, transform_cvar::Float64, scaling::Float64)
    inv = cost_df[!,"Investment_Cost"] ./ scaling
    cvar = cost_df[!,"CVaR_OpCost"] ./ scaling
    ev = cost_df[!,"EV_OpCost"] ./ scaling

    # Drop solutions the Benders mapping recording pass let through despite violating the
    # Transformed budget constraint (it only checks the Benders convergence gap, not the
    # budget itself) - see algorithm.jl:432. Filtered against the loosest (untightened)
    # budget since results from every tightening level are combined here.
    transformed_cost = transform_sys .* (inv .+ ev) .+ transform_cvar .* (cvar)
    cutoff = budget + settings["Mapping Gap Threshold"]
    @info("Budget in sampling: $(budget); cutoff: $(cutoff)")
    keep = transformed_cost .<= cutoff

    n_dropped = length(keep) - sum(keep)
    if n_dropped > 0
        @info("Dropping $(n_dropped) of $(length(keep)) mapped solutions with transformed cost above budget + Mapping Gap Threshold ($(cutoff)) before Delaunay sampling.")
    end



    points = points[keep, :]
    inv, cvar, ev = inv[keep], cvar[keep], ev[keep]
    cost_points = Matrix([inv+ev inv+cvar])

    return sample_interior_simplex(points, cost_points, num_samples, settings)
end
