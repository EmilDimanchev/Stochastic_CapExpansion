###########

# Interpolation evaluation functions

#######

function evaluate_interpolates(SPs::Array{Model, 3}, inputs::Dict, settings::Dict, interpolate_df::DataFrame, results_folder::String)
    # result containers
    results_cap = []
    results_syscost = []
    results_emissions = []

    P_s = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]
    VaR_percent = settings["Value-at-Risk percent"]
    Z = inputs["Number of zones"]
    costs_by_resource = inputs["Investment costs"]
    cost_of_lines = inputs["CAPEX per MW"]
    n_lines = length(cost_of_lines)

    # isolate capacity vectors
    index_indx = findfirst(==("Index"), names(interpolate_df))
    labels = Vector(interpolate_df[:, index_indx])
    select!(interpolate_df, Not(:Index))
    investment_indx = findfirst(==("Investment_Cost"), names(interpolate_df))
    cap_vectors = interpolate_df[:, 1:investment_indx-n_lines-1]
    line_vectors = interpolate_df[:, investment_indx-n_lines:investment_indx-1]
    # get costs
    cap_inv_costs = cap_vectors .* reshape(costs_by_resource, 1, :)
    line_inv_costs = line_vectors .* reshape(cost_of_lines, 1, :)
    inv_costs = hcat(cap_inv_costs, line_inv_costs)
    tot_inv_costs = sum.(eachrow(Matrix(inv_costs)))
    
    col_by_zone = collect(findall(x -> occursin("z"*string(z), x), names(inv_costs)) for z in 1:Z)
    cost_by_zone = [sum.(eachrow(inv_costs[:, col_by_zone[z]])) for z in 1:Z]

    for i in 1:nrow(interpolate_df)
        @info("Evaluating interpolate ", labels[i])
        caps = Vector(cap_vectors[i, :])
        lines = Vector(line_vectors[i, :])
        outputs_sp = run_all_subproblems(SPs, inputs, settings, caps, lines; minimal_payload=false)
        ev, cvar = evaluate_subproblems(outputs_sp, P_s, P_f, P_k, VaR_percent)

        costs_by_zone_it = [cost_by_zone[z][i] for z in 1:Z]

        df_cap, df_syscost, df_emissions = make_results_mapping_dfs(caps, lines, outputs_sp, tot_inv_costs[i], costs_by_zone_it, cvar, ev, inputs, settings)
        push!(results_cap, df_cap)
        push!(results_syscost, df_syscost)
        push!(results_emissions, df_emissions)
    end
    write_exploration_results!(results_cap, results_syscost, results_emissions, results_folder, labels, "Interpolates"; mapping=true)
end

function evaluate_subproblems(outputs_sp::Array{Dict{String, Any}, 3}, P_s::Vector, P_f::Vector, P_k::Vector, VaR_percent::Float64)
    S = length(P_s)
    F = length(P_f)
    K = length(P_k)

    sp_obj_per_iter = reshape([outputs_sp[s,f,k]["SP objective"] for s in 1:S, f in 1:F, k in 1:K], (S, F, K))
    ev = sum(P_s[s]*P_f[f]*P_k[k]*sp_obj_per_iter[s,f,k] for s in 1:S, f in 1:F, k in 1:K)
    cvar = compute_cvar(sp_obj_per_iter, P_s, P_f, P_k, VaR_percent)
    return ev, cvar
end

#==========================================================================================

Planner design CSVs (e.g. experiments/designs/run13_planner_designs.csv)

One design per row, keyed by point_id, capacities as <Tech>_Z<z> columns in GW plus a
Transmission column, alongside the planner's predicted metrics. read_design_csv maps the
capacity columns by name onto capacity_variable_names(inputs); every other column is kept as
Design_<name> metadata. evaluate_design_csv runs the designs through run_subproblems_batch and
builds Summary rows exactly as evaluate_picks does, with the metadata appended.

==========================================================================================#

const DESIGN_CSV_TECHS = Dict("GasCC" => "Gas", "GasCT" => "Gas turbine", "Wind" => "Wind", "Solar" => "Solar",
                              "Nuclear" => "Nuclear", "Coal" => "Coal", "Battery" => "Battery")

# Design CSV column for each capacity variable, in capacity_variable_names(inputs) order
function design_csv_columns(inputs::Dict, csv_names::AbstractVector)
    by_resource = Dict{String,String}()
    for col in csv_names
        m = match(r"^([A-Za-z0-9]+)_Z(\d+)$", col)
        (m === nothing || !haskey(DESIGN_CSV_TECHS, m[1])) && continue
        resource = DESIGN_CSV_TECHS[m[1]] * "-z" * m[2]
        resource in inputs["Resources"] || error("Design CSV column $col maps to $resource, which is not in inputs[\"Resources\"].")
        by_resource[resource] = col
    end
    missing_res = setdiff(inputs["Resources"], keys(by_resource))
    isempty(missing_res) || error("Design CSV has no capacity column for $(missing_res).")
    L = inputs["Number of lines"]
    line_cols = L == 1 ? ["Transmission"] : ["Transmission_" * string(l) for l in 1:L]
    absent = setdiff(line_cols, csv_names)
    isempty(absent) || error("Design CSV is missing line column(s) $(absent).")
    return vcat([by_resource[r] for r in inputs["Resources"]], line_cols)
end

# (labels, X, meta): point_ids, capacity vectors in MW (resources then lines), and the
# remaining columns renamed Design_<name>, one row per design
function read_design_csv(df::DataFrame, inputs::Dict; unit_scale::Real = 1000.0, inv_rel_tol::Real = 1e-3)
    "point_id" in names(df) || error("Design CSV has no point_id column.")
    cols = design_csv_columns(inputs, names(df))
    labels = String[string(v) for v in df.point_id]
    allunique(labels) || error("Design CSV has duplicate point_ids.")
    R = length(inputs["Resources"])
    ub = vcat(inputs["Capacity upper bounds"], inputs["Max expansion"])
    X = Vector{Vector{Float64}}(undef, nrow(df))
    for i in 1:nrow(df)
        x = Float64[df[i, c] for c in cols] .* unit_scale
        neg = findall(<(0.0), x)
        for j in neg
            x[j] >= -1e-6 || error("Design $(labels[i]): $(cols[j]) = $(x[j]) MW is negative.")
            @warn("Design $(labels[i]): clipping $(cols[j]) = $(x[j]) MW to 0.")
            x[j] = 0.0
        end
        over = findall(x .> ub)
        isempty(over) || error("Design $(labels[i]): $(cols[over]) exceed their upper bounds $(ub[over]) MW.")
        if "Investment_B" in names(df)
            inv, _ = investment_cost(x, inputs)
            claimed = df[i, "Investment_B"] * 1e9
            abs(inv - claimed) <= inv_rel_tol * abs(claimed) ||
                @warn("Design $(labels[i]): investment cost $(inv) differs from Investment_B*1e9 = $(claimed); check units/columns.")
        end
        X[i] = x
    end
    meta = select(df, Not(vcat(cols, "point_id")))
    rename!(n -> "Design_" * n, meta)
    @info("Read $(length(X)) designs ($(R) resources + $(length(cols) - R) lines, ×$(unit_scale) to MW).")
    return labels, X, meta
end

read_design_csv(path::String, inputs::Dict; kwargs...) = read_design_csv(CSV.read(path, DataFrame), inputs; kwargs...)

# Evaluates the designs (all, or the row numbers in `rows`) in batches of batch_size, writing
# results_folder/Summary_<summary_name>.csv after every batch as a checkpoint.
function evaluate_design_csv(SPs, inputs::Dict, settings::Dict, design_file::String, results_folder::String;
                             batch_size::Int = 25, unit_scale::Real = 1000.0, rows = nothing,
                             summary_name::String = splitext(basename(design_file))[1])
    labels, X, meta = read_design_csv(design_file, inputs; unit_scale = unit_scale)
    idx = rows === nothing ? collect(eachindex(X)) : collect(rows)
    R = length(inputs["Resources"])
    p = scenario_probabilities(inputs)
    α = settings["Value-at-Risk percent"]
    rows_cap, rows_syscost, rows_emissions = DataFrame[], DataFrame[], DataFrame[]
    done = String[]
    n_batches = cld(length(idx), batch_size)
    for (b, batch) in enumerate(Iterators.partition(idx, batch_size))
        t = @elapsed sps = run_subproblems_batch(SPs, inputs, settings, [X[i][1:R] for i in batch], [X[i][R+1:end] for i in batch]; minimal_payload = false)
        for (i, sp) in zip(batch, sps)
            inv, inv_by_zone = investment_cost(X[i], inputs)
            fs = final_solve_from_sp_outputs(labels[i], X[i], sp, p, inv)
            ev, cvar = solve_ev(fs), solve_cvar(fs, α)
            df_cap, df_syscost, df_emissions = make_results_mapping_dfs(X[i][1:R], X[i][R+1:end], sp, inv, inv_by_zone, cvar, ev, inputs, settings)
            df_syscost = hcat(df_syscost, meta[i:i, :])
            df_syscost[!, :Design_Solve_OK] = [fs.ok]
            push!(rows_cap, df_cap); push!(rows_syscost, df_syscost); push!(rows_emissions, df_emissions)
            push!(done, labels[i])
        end
        write_exploration_results!(rows_cap, rows_syscost, rows_emissions, results_folder, done, summary_name; mapping = true)
        @info("Design batch $(b)/$(n_batches): $(length(batch)) designs in $(round(t; digits = 1)) s ($(length(done))/$(length(idx)) done).")
    end
    return done
end
