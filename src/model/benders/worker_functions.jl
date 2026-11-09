export build_subproblems, set_capacity_parameters!, run_subproblem
using Distributed, DistributedArrays, JuMP, Gurobi

#==========

Subproblems are stored as an S × F × K × B array: one model per demand, gas and weather scenario
and per time block b. "Subproblems per scenario" in settings sets B, and each block is an equal
run of whole representative periods (B = 1 is one subproblem per scenario). No constraint links
two periods (storage wraps and the energy limit apply within each period), so splitting a
scenario into blocks is exact: its objective and duals are the sums over its blocks
(merge_block_results). Everything downstream
of run_all_subproblems sees one result per scenario, as before.

==========#

# (s, f, k, b): demand, gas and weather scenario, then time block
const SubproblemKey = NTuple{4, Int}
const WORKER_SP_CACHE = Dict{SubproblemKey, Model}()
const WORKER_DATA_CACHE = Dict{Symbol, Any}()
const MASTER_WORKER_SP_KEYS = Dict{Int, Set{SubproblemKey}}()
const MASTER_SCENARIO_TO_WORKER = Dict{SubproblemKey, Int}()
const WORKER_SOLVE_COUNTER = Ref(0)

function _worker_memory_debug_enabled()
    settings = get(WORKER_DATA_CACHE, :settings, nothing)
    return settings isa Dict && get(settings, "Memory debug flag", false)
end

function _worker_memory_log_frequency()
    settings = get(WORKER_DATA_CACHE, :settings, nothing)
    if !(settings isa Dict)
        return 1
    end
    raw = get(settings, "Worker memory log frequency", 1)
    return max(1, Int(raw))
end

function _worker_log_memory!(label::String, scenario_key::SubproblemKey; payload_size_mb::Union{Nothing, Float64}=nothing)
    if !_worker_memory_debug_enabled()
        return nothing
    end
    free_mem_mb = round(Sys.free_memory() / 1024^2; digits = 2)
    cache_size = length(WORKER_SP_CACHE)
    if isnothing(payload_size_mb)
        @info("[worker $(myid())] $(label) scenario=$(scenario_key) cache=$(cache_size) free=$(free_mem_mb) MiB")
    else
        @info("[worker $(myid())] $(label) scenario=$(scenario_key) cache=$(cache_size) free=$(free_mem_mb) MiB payload=$(payload_size_mb) MiB")
    end
    return nothing
end

function _set_worker_problem_data!(inputs, settings)
    WORKER_DATA_CACHE[:inputs] = inputs
    WORKER_DATA_CACHE[:settings] = settings
    return nothing
end

function _clear_worker_subproblem_cache!()
    empty!(WORKER_SP_CACHE)
    GC.gc(false)
    return nothing
end

function _build_and_cache_subproblem!(scenario_key::SubproblemKey)
    inputs, settings = WORKER_DATA_CACHE[:inputs], WORKER_DATA_CACHE[:settings]
    s, f, k, b = scenario_key
    SP = build_subproblem(inputs, settings, [s, f, k]; hours = subproblem_blocks(inputs, settings)[b])
    WORKER_SP_CACHE[scenario_key] = SP
    return nothing
end

function _run_cached_subproblem!(scenario_key::SubproblemKey, capacity::Vector{Float64}, capacity_line::Vector{Float64}, minimal_payload::Bool)
    SP = WORKER_SP_CACHE[scenario_key]
    set_parameter_value.(SP[:x], capacity)
    set_parameter_value.(SP[:x_line], capacity_line)
    result = run_subproblem(SP, WORKER_DATA_CACHE[:inputs], WORKER_DATA_CACHE[:settings]; minimal_payload=minimal_payload)
    return result
end

# Several capacity points on this worker's own scenarios, one after the other. Each cached
# model only ever sees one point at a time and warm-starts from the previous point's basis.
# Returns results[point][key index].
function _run_cached_batch!(scenario_keys::Vector{SubproblemKey}, capacities::Vector{Vector{Float64}}, capacity_lines::Vector{Vector{Float64}}, minimal_payload::Bool)
    return [[_run_cached_subproblem!(key, cap, line, minimal_payload) for key in scenario_keys] for (cap, line) in zip(capacities, capacity_lines)]
end

# Hour ranges of the time blocks each scenario is split into. "Subproblems per scenario" = N gives
# N blocks of equal numbers of whole representative periods (missing = 1, one block of every period).
# N must divide the period count, so a block never cuts through a period and blocks are equally sized.
function subproblem_blocks(inputs, settings)
    first_hours, last_hours = inputs["Set of first periods"], inputs["Set of last periods"]
    n_periods = length(first_hours)
    valid = filter(d -> n_periods % d == 0, 1:n_periods)
    raw = get(settings, "Subproblems per scenario", 1)
    (raw isa Integer && raw in valid) || error("Subproblems per scenario is $(repr(raw)), but each scenario has $(n_periods) representative periods, " *
        "which must split evenly into whole periods per subproblem. Valid values: $(join(valid, ", ")).")
    per_block = n_periods ÷ raw
    return [first_hours[first(p)]:last_hours[last(p)] for p in Iterators.partition(1:n_periods, per_block)]
end

# Demand, gas and weather scenario counts (S, F, K) of a subproblem array, without its time blocks
scenario_dims(SPs::Array{Model,4}) = size(SPs)[1:3]

function build_all_subproblems(inputs, settings)
    S = inputs["Number of demand scenarios"]
    F = inputs["Number of gas price scenarios"]
    K = inputs["Number of weather scenarios"]
    blocks = subproblem_blocks(inputs, settings)
    B = length(blocks)
    SPs = Array{Model,4}(undef, S, F, K, B)
    # A single scenario stays on the main process even when it has several blocks: a parallel build
    # clears every worker's cache, which would drop a stochastic set already living there
    # (run_base_mga builds its one-scenario set after the full one).
    if settings["Parallel flag"] && nworkers() > 0 && (S > 1 || F > 1 || K > 1)
        @info "Building $(S*F*K) scenarios × $(B) time blocks in parallel using $(nworkers()) distributed workers..."
        pids = workers()
        empty!(MASTER_SCENARIO_TO_WORKER)
        @sync for pid in pids
            @async begin
                remotecall_wait(_clear_worker_subproblem_cache!, pid)
                remotecall_wait(_set_worker_problem_data!, pid, inputs, settings)
                worker_set = get!(MASTER_WORKER_SP_KEYS, pid, Set{SubproblemKey}())
                empty!(worker_set)
            end
        end
        # Round-robin over (scenario, block), so one scenario's blocks are spread across workers
        scenario_keys = [(s, f, k, b) for s in 1:S, f in 1:F, k in 1:K, b in 1:B]
        @sync for (idx, scenario_key) in enumerate(scenario_keys)
            @async begin
                s, f, k, b = scenario_key
                SPs[s,f,k,b] = build_subproblem(inputs, settings, [s, f, k]; hours = blocks[b])
                pid = pids[mod1(idx, length(pids))]
                remotecall_wait(_build_and_cache_subproblem!, pid, scenario_key)
                push!(MASTER_WORKER_SP_KEYS[pid], scenario_key)
                MASTER_SCENARIO_TO_WORKER[scenario_key] = pid
            end
        end
    else
        if settings["Parallel flag"] && nworkers() == 0
            @warn("Parallel flag is true but no distributed workers are available. Falling back to serial subproblem build.")
        elseif (S == 1 && F == 1 && K == 1)
            @info("One Scenario selected, building its $(B) time block(s) on main process")
        end
        for s in 1:S, f in 1:F, k in 1:K, b in 1:B
            SPs[s,f,k,b] = build_subproblem(inputs, settings, [s,f,k]; hours = blocks[b])
        end
    end
    return SPs
end

# The SP works on hours (a whole number of representative periods, see subproblem_blocks), indexed
# locally as 1:length(hours). Each period keeps its own storage wrap and energy limit.
function build_subproblem(inputs, settings, scenario_index::AbstractVector; hours::UnitRange{Int} = 1:inputs["Number of periods"])
    # Set Scenario index
    s = scenario_index[1]
    f = scenario_index[2]
    k = scenario_index[3]

    # Periods inside the block, as local hour indices. A block that cuts through a period would
    # break that period's storage wrap and energy limit, and the split would no longer be exact.
    all_first, all_last = inputs["Set of first periods"], inputs["Set of last periods"]
    (first(hours) in all_first && last(hours) in all_last) || error("Subproblem hours $(hours) don't start and end on representative period boundaries.")
    block_periods = findall(i -> all_first[i] >= first(hours) && all_last[i] <= last(hours), eachindex(all_first))
    offset = first(hours) - 1

    demand_scenario = inputs["Demand"][hours,:,k]
    variable_costs_scenario = inputs["Variable costs"][:,f]
    availability_scenario = inputs["Generation availability"][hours,:,k]
    period_weights_scenario = inputs["Period weights"][hours,k]
    demand_adder_scenario = inputs["Demand adders"][s,:]

     # Numerical settings
    scaling_factor_demand = settings["Scaling factor demand"]
    scaling_factor_cost = settings["Scaling factor cost"] 

    # Model
    ED = Model()
    ED.ext[:Variable_costs_scenario] = variable_costs_scenario
    ED.ext[:Availability_scenario] = availability_scenario
    ED.ext[:Period_weights_scenario] = period_weights_scenario
    ED.ext[:Demand_adder_scenario] = demand_adder_scenario
    if settings["Solver"] == "HiGHS"
        set_optimizer(ED, HiGHS.Optimizer)
        set_optimizer_attribute(ED, "primal_feasibility_tolerance", 1e-3)
        set_optimizer_attribute(ED, "optimality_tolerance", 1e-3)
        # Force crossover on so a vertex/basis is always produced: x/x_line are Parameter
        # variables (see below), so re-solves after set_parameter_value! are pure bound
        # changes on an otherwise-static LP - a basis lets HiGHS warm-start the next
        # solve instead of re-deriving a starting point from scratch every iteration.
        set_optimizer_attribute(ED, "run_crossover", "on")
    elseif settings["Solver"] == "Gurobi"
        set_optimizer(ED, Gurobi.Optimizer)
        set_optimizer_attribute(ED, "OutputFlag", 0)
        set_optimizer_attribute(ED, "Crossover", 1)
        # Automatic (-1) rather than hard-pinned barrier: with Crossover=1 a basis is
        # always available from the previous solve, and Gurobi's automatic method
        # selection favors warm-started simplex continuation when one exists while still
        # falling back to barrier/concurrent for a cold or very large LP - this adapts
        # across problem sizes rather than committing to one method.
        set_optimizer_attribute(ED, "Method", -1)
        set_optimizer_attribute(ED, "BarConvTol", 1e-5)
        set_optimizer_attribute(ED, "OptimalityTol", 1e-5)
        set_optimizer_attribute(ED, "FeasibilityTol", 1e-5)
    end
    set_silent(ED)

     # ~~~~
    # Load inputs
    # ~~~

    # Settings
    storage_flag = inputs["Storage flag"]
    crm_flag = settings["CRM flag"]

    # Parameters
    availability = availability_scenario
    co2_factors = inputs["CO2 emission intensities"]

    # Storage specific
    if storage_flag
        N = inputs["Storage power to energy ratio"]
        # Efficiency
        F_dch = inputs["Storage discharge loss"]
        F_ch = inputs["Storage charge loss"]

        first_periods = all_first[block_periods] .- offset
        last_periods = all_last[block_periods] .- offset
    end

    # Representative periods in this block
    n_periods = length(block_periods)
    t_weights = period_weights_scenario
    zone_map = inputs["Zone map"]
    res_zone_map = inputs["Resource zone map"]
    line_capacity_existing = inputs["Line capacities"] # MW, in place before any expansion

    # ~~~
    # Build model
    # ~~~

    # Sets
    T = length(hours)
    G = inputs["Number of generation resources"]
    O =  inputs["Number of storage resources"]
    storage_flag = O > 0
    R = G + O
    L = inputs["Number of lines"]
    Z = inputs["Number of zones"]
    # Only each zone's own resources and lines. Summing over all of them with the maps as
    # coefficients stores a zero term per pair, which gets expensive with many zones.
    zone_generators = [findall(>(0), res_zone_map[1:G, z]) for z in 1:Z]
    zone_storage = [findall(>(0), res_zone_map[G+1:R, z]) for z in 1:Z]
    zone_lines = [findall(!=(0), zone_map[:, z]) for z in 1:Z]

    # ~~~
    # Parameters
    # ~~~
    
    # Scale parameters
    
    demand_base = demand_scenario./scaling_factor_demand
    demand_add = demand_adder_scenario./scaling_factor_demand
    cost_var = variable_costs_scenario./scaling_factor_cost
    # Demand
    @expression(ED, demand[t in 1:T, z in 1:Z], demand_base[t, z] + demand_add[z])
    max_nse = settings["Max nse"] # percentage of demand that can be curtailed in each segments
    nse_segs = length(max_nse)
    price_nse = settings["Cost of nonserved energy"]./scaling_factor_cost

    # ~~~
    # Model formulation
    # ~~~

    # Generation
    @variable(ED, g[r in 1:G, t in 1:T] >= 0) # MWh
    @variable(ED, x[r in 1:R] in Parameter(0.0)) # Capacity, MW
    @variable(ED, x_line[l in 1:L] in Parameter(0.0)) # Line capacity, MW
    # Non-servED energy
    @variable(ED, y[seg in 1:nse_segs, t in 1:T, z in 1:Z] >= 0) # $/MWh
    @expression(ED, cost_nse[t in 1:T, z in 1:Z], sum(price_nse[seg]*y[seg,t,z] for seg in 1:nse_segs))

    if storage_flag
        # One set per storage resource (rows of Resources_storage.csv)
        @variable(ED, z_dch[o in 1:O, t in 1:T] >= 0)
        @variable(ED, z_ch[o in 1:O, t in 1:T] >= 0)
        # State of charge
        @variable(ED, e[o in 1:O, t in 1:T] >= 0)
        # Primal feasibility

        for o in 1:O, i in 1:n_periods
            first_t, last_t = first_periods[i], last_periods[i]
            # Wrap first and last periods - start
            @constraint(ED, e[o, first_t] == e[o, last_t] - (1/F_dch[o])*z_dch[o, first_t] + F_ch[o]*z_ch[o, first_t])
            @constraint(ED, z_dch[o, first_t] <= e[o, last_t])
            # Wrap first and last periods - end
            @constraint(ED, [t in first_t+1:last_t], e[o, t] == e[o, t-1] - (1/F_dch[o])*z_dch[o, t] + F_ch[o]*z_ch[o, t])
            @constraint(ED, [t in first_t+1:last_t], z_dch[o, t] <= e[o, t-1])
        end


        # Energy balance for the remaining periods
        @constraint(ED, energy_limit[o in 1:O, t in 1:T], e[o, t] <= (1/N[o])*x[G+o])
        @constraint(ED, charge_limit_rate[o in 1:O, t in 1:T], z_ch[o, t] <= x[G+o])

        @constraint(ED, discharge_limit_rate[o in 1:O, t in 1:T], z_dch[o, t] <= x[G+o])
        @constraint(ED, charge_discharge_balance[o in 1:O, t in 1:T], z_dch[o, t] + z_ch[o, t] <= x[G+o])
    end

    # Load shedding
    # Maximum non served energy (primal feasibility)
    @constraint(ED, max_load_shedding[seg in 1:nse_segs,t in 1:T, z in 1:Z], max_nse[seg]*demand[t,z] - y[seg,t,z] >= 0)
    @expression(ED, total_nse[t in 1:T, z in 1:Z], sum(y[seg,t,z] for seg in 1:nse_segs))
    # Consumer surplus including demand shock for reporting purposes
    @expression(ED, served_demand[seg in 1:nse_segs, t in 1:T, z in 1:Z], max_nse[seg]*demand[t,z] - y[seg,t,z])
    @expression(ED, total_benefit[t in 1:T, z in 1:Z], price_nse[1]*demand_base[t,z] + sum(price_nse[seg]*served_demand[seg,t,z] for seg in 2:nse_segs))

    # Capacity limit on generation
    @constraint(ED, capacity_limit[r in 1:G, t in 1:T], g[r,t] <= x[r]*availability[t,r])
    # Energy-limited generators (e.g. hydro): output over each representative period can't exceed
    # the fraction of what full capacity would produce, however it's shaped within the period
    energy_fractions = inputs["Energy limit fractions"]
    period_hours = [(all_first[i] - offset):(all_last[i] - offset) for i in block_periods]
    @constraint(ED, generation_energy_limit[r in findall(<(1.0), energy_fractions), i in 1:n_periods],
        sum(g[r,t] for t in period_hours[i]) <= energy_fractions[r]*length(period_hours[i])*x[r])
    if Z > 1
        @variable(ED, flow[t in 1:T, l in 1:L])
        @constraint(ED, transmission_limit[t in 1:T, l in 1:L], flow[t,l] <= (line_capacity_existing[l] + x_line[l])/scaling_factor_demand)
        @constraint(ED, transmission_limit_negative[t in 1:T, l in 1:L], flow[t,l] >= -(line_capacity_existing[l] + x_line[l])/scaling_factor_demand)
    end

    # Power balance constraint
    @expression(ED, supply[t in 1:T, z in 1:Z], AffExpr(0.0))
    for t in 1:T, z in 1:Z
        for r in zone_generators[z]
            add_to_expression!(supply[t,z], res_zone_map[r,z], g[r,t])
        end
        if storage_flag
            for o in zone_storage[z]
                add_to_expression!(supply[t,z], res_zone_map[G+o,z], z_dch[o,t])
                add_to_expression!(supply[t,z], -res_zone_map[G+o,z], z_ch[o,t])
            end
        end
        if Z > 1
            for l in zone_lines[z]
                add_to_expression!(supply[t,z], zone_map[l,z], flow[t,l])
            end
        end
    end
    @constraint(ED, power_balance[t in 1:T, z in 1:Z], supply[t,z] + total_nse[t,z] - demand[t,z] == 0)

    # Hourly Zonal CRM
    if crm_flag
        crm_margin = inputs["CRM Margin"]
        crm_price_cap = inputs["CRM Price Cap"]./scaling_factor_cost
        crm_derate = inputs["CRM Derating Factors"]
        gen_resource_names = inputs["Generation resources"]
        derate_gen = [crm_derate[gen_resource_names[r]] for r in 1:G]

        # Eligible capacity: derated generation availability, derated storage net discharge,
        # derated import/export flows, and voluntary curtailment (segments 2+, excludes emergency shedding)
        @expression(ED, eCapResMargin[t in 1:T, z in 1:Z],
            sum(derate_gen[r]*x[r]*availability[t,r]*res_zone_map[r,z] for r in zone_generators[z]; init = AffExpr(0.0))
            + sum(y[seg,t,z] for seg in 2:nse_segs))
        for t in 1:T, z in 1:Z
            if storage_flag
                for o in zone_storage[z]
                    add_to_expression!(eCapResMargin[t,z], res_zone_map[G+o,z], z_dch[o,t])
                    add_to_expression!(eCapResMargin[t,z], -res_zone_map[G+o,z], z_ch[o,t])
                end
            end
            if Z > 1
                for l in zone_lines[z]
                    add_to_expression!(eCapResMargin[t,z], zone_map[l,z], flow[t,l])
                end
            end
        end

        @variable(ED, vCRMSlack[t in 1:T, z in 1:Z] >= 0)
        @constraint(ED, cCapacityResMargin[t in 1:T, z in 1:Z],
            eCapResMargin[t,z] + vCRMSlack[t,z] >= (1 + crm_margin[z])*demand[t,z])
        @expression(ED, crm_slack_cost_by_zone[z in 1:Z], sum(t_weights[t]*crm_price_cap[z]*vCRMSlack[t,z] for t in 1:T))
    else
        @expression(ED, crm_slack_cost_by_zone[z in 1:Z], AffExpr(0.0))
    end
    @expression(ED, crm_slack_cost, sum(crm_slack_cost_by_zone[z] for z in 1:Z))


    # Weighted the same way as the objective, so the zones sum to the SP objective
    @expression(ED, cost_by_zone[z in 1:Z], sum(t_weights[t]*g[r,t]*res_zone_map[r,z]*cost_var[r] for r in zone_generators[z], t in 1:T; init = AffExpr(0.0)) + sum(t_weights[t]*cost_nse[t,z] for t in 1:T) + crm_slack_cost_by_zone[z])
    @expression(ED, emissions_by_zone[z in 1:Z], sum(t_weights[t]*g[r,t]*res_zone_map[r,z]*co2_factors[r] for r in zone_generators[z], t in 1:T; init = AffExpr(0.0)))
    # ~~~
    # Objective function
    # ~~~

    @objective(ED, Min, sum(t_weights[t]*g[r,t]*cost_var[r] for r in 1:G, t in 1:T) + sum(t_weights[t]*cost_nse[t,z] for t in 1:T, z in 1:Z) + crm_slack_cost)

    return ED
end


function set_capacity_parameters!(SPs::Array{Model,4}, capacity::Vector{Float64}, capacity_line::Vector{Float64})
    for SP in SPs
        set_parameter_value.(SP[:x], capacity)
        set_parameter_value.(SP[:x_line], capacity_line)
    end
end

function run_subproblem(ED::Model, inputs, settings; minimal_payload::Bool=false)
    output = Dict{String, Any}()

    optimize!(ED)
    _log_solver_work!("SP", ED, settings)
    if termination_status(ED) != MOI.OPTIMAL
        @warn("Subproblem did not solve to optimality. Status: $(termination_status(ED))")
        output["coeff"] = 0
    else
        output["coeff"] = 1
    end
    
    output_full = write_subproblem_results(ED,inputs, settings; minimal_payload=minimal_payload)
    output = merge(output, output_full)

    return output
end

# How each write_subproblem_results key combines across a scenario's time blocks (in block order).
# Totals add up, the solve status is the worst block's, and hourly series are stacked along time.
const BLOCK_SUM_KEYS = ("SP objective", "capacity dual", "line dual", "Cost by zone", "Emissions by zone", "Emissions")
const BLOCK_TIME_ROW_KEYS = ("Power price", "Load shedding", "Net Flow", "CRM dual", "CRM slack")        # T × ·
const BLOCK_TIME_COLUMN_KEYS = ("Generation", "Storage charging", "Storage discharging")                # · × T

# One scenario's result from the results of its blocks. Every key needs a rule above, so a new
# subproblem output can't be silently dropped.
function merge_block_results(outs::AbstractVector)
    length(outs) == 1 && return outs[1]
    merged = Dict{String, Any}()
    for key in keys(outs[1])
        parts = [o[key] for o in outs]
        if key == "coeff"
            merged[key] = minimum(parts)
        elseif key in BLOCK_SUM_KEYS
            merged[key] = sum(parts)
        elseif key in BLOCK_TIME_ROW_KEYS
            merged[key] = reduce(vcat, parts)
        elseif key in BLOCK_TIME_COLUMN_KEYS
            merged[key] = reduce(hcat, parts)
        else
            error("merge_block_results: no rule for subproblem output \"$(key)\". Add it to one of the BLOCK_*_KEYS lists.")
        end
    end
    return merged
end

# Per-block results (S × F × K × B) to per-scenario results (S × F × K)
function merge_block_results(block_results::Array{<:Any,4})
    S, F, K, _ = size(block_results)
    return [merge_block_results(view(block_results, s, f, k, :)) for s in 1:S, f in 1:F, k in 1:K]
end

_run_in_parallel(SPs::Array{Model,4}, settings) = settings["Parallel flag"] && nworkers() > 0 && scenario_dims(SPs) != (1, 1, 1)

function run_all_subproblems(SPs::Array{Model,4}, inputs, settings, capacity::Vector{Float64}, capacity_line::Vector{Float64}; minimal_payload::Bool=get(settings, "Minimal worker payload flag", true))

    block_results = Array{Dict{String, Any},4}(undef, size(SPs)...)
    if _run_in_parallel(SPs, settings)
        pids = workers()
        scenario_keys = [Tuple(I) for I in CartesianIndices(SPs)]

        @sync for scenario_key in scenario_keys
            @async begin
                if !haskey(MASTER_SCENARIO_TO_WORKER, scenario_key)
                    error("Subproblem $(scenario_key) has no worker assignment. Call build_all_subproblems() before run_all_subproblems().")
                end
                pid = MASTER_SCENARIO_TO_WORKER[scenario_key]
                if !(pid in pids)
                    error("Assigned worker $(pid) for subproblem $(scenario_key) is not active. Rebuild subproblems with build_all_subproblems().")
                end
                worker_keys = get(MASTER_WORKER_SP_KEYS, pid, Set{SubproblemKey}())
                if !(scenario_key in worker_keys)
                    error("Subproblem $(scenario_key) is not cached on worker $(pid). Call build_all_subproblems() before run_all_subproblems().")
                end
                block_results[scenario_key...] = @fetchfrom pid _run_cached_subproblem!(scenario_key, capacity, capacity_line, minimal_payload)

            end
        end
        @everywhere GC.gc()
    else
        if settings["Parallel flag"] && nworkers() == 0
            @warn("Parallel flag is true but no distributed workers are available. Falling back to serial subproblem solves.")
        elseif scenario_dims(SPs) == (1, 1, 1)
            @info("One scenario loaded, running on main process")
        end
        set_capacity_parameters!(SPs, capacity, capacity_line)
        for I in CartesianIndices(SPs)
            block_results[I] = run_subproblem(SPs[I], inputs, settings; minimal_payload=minimal_payload)
        end
    end
    return merge_block_results(block_results)
end

# run_all_subproblems for a batch of capacity points. Subproblems stay on the workers that
# own them (build_all_subproblems); each worker gets ONE call that loops over the batch for its
# own subproblems, so no two in-flight calls share a cached model and no per-point GC/dispatch.
# Returns one S × F × K results array per point, identical to run_all_subproblems' output.
function run_subproblems_batch(SPs::Array{Model,4}, inputs, settings, capacities::Vector{Vector{Float64}}, capacity_lines::Vector{Vector{Float64}}; minimal_payload::Bool=get(settings, "Minimal worker payload flag", true))
    n = length(capacities)
    length(capacity_lines) == n || error("run_subproblems_batch: $(n) capacity vectors but $(length(capacity_lines)) line vectors.")
    if !_run_in_parallel(SPs, settings)
        return [run_all_subproblems(SPs, inputs, settings, capacities[i], capacity_lines[i]; minimal_payload=minimal_payload) for i in 1:n]
    end
    block_results = [Array{Dict{String, Any},4}(undef, size(SPs)...) for _ in 1:n]
    owned = Dict(pid => sort!(collect(keys_set)) for (pid, keys_set) in MASTER_WORKER_SP_KEYS if pid in workers() && !isempty(keys_set))
    n_owned = sum(length, values(owned); init = 0)
    n_owned == length(SPs) || error("run_subproblems_batch: workers own $(n_owned) of $(length(SPs)) subproblems. Call build_all_subproblems() first.")
    @sync for (pid, scenario_keys) in owned
        @async begin
            out = remotecall_fetch(_run_cached_batch!, pid, scenario_keys, capacities, capacity_lines, minimal_payload)
            for i in 1:n, (j, key) in enumerate(scenario_keys)
                block_results[i][key...] = out[i][j]
            end
        end
    end
    @everywhere GC.gc()
    return [merge_block_results(r) for r in block_results]
end

function write_subproblem_results(ED::Model, inputs, settings; minimal_payload::Bool=false)
    output = Dict{String, Any}()

    
    availability_scenario = ED.ext[:Availability_scenario]
    period_weights_scenario = ED.ext[:Period_weights_scenario]

    # Settings
    storage_flag = inputs["Storage flag"]
    crm_flag = settings["CRM flag"]

    # Parameters
    co2_factors = inputs["CO2 emission intensities"]
    
    # Representative periods
    t_weights = period_weights_scenario


    # Sets (T is this block's hours, see build_subproblem)
    T = length(t_weights)
    G = inputs["Number of generation resources"]
    O =  inputs["Number of storage resources"]
    storage_flag = O > 0
    L = inputs["Number of lines"]
    Z = inputs["Number of zones"]


    max_nse = settings["Max nse"] # percentage of demand that can be curtailed in each segments
    nse_segs = length(max_nse)
    
    
    output["capacity dual"] = dual.(ParameterRef.(ED[:x])).*settings["Scaling factor cost"]
    output["line dual"] = dual.(ParameterRef.(ED[:x_line])).*settings["Scaling factor cost"]

    output["SP objective"] = objective_value(ED).*(settings["Scaling factor cost"]*settings["Scaling factor demand"])
    if minimal_payload
        return output
    end
    output["Cost by zone"] = value.(ED[:cost_by_zone]).*settings["Scaling factor cost"]
    output["Emissions by zone"] = value.(ED[:emissions_by_zone])
    output["Power price"] = dual.(ED[:power_balance]).*settings["Scaling factor cost"]./t_weights
    nse_dual = dual.(ED[:max_load_shedding]).*settings["Scaling factor cost"]
    #output["Max NSE dual"] = sum(max_nse[seg]*ED[:demand][t]*nse_dual[seg,t] for seg in 1:nse_segs, t in 1:T)
    generation = value.(ED[:g])
    output["Emissions"] = sum(t_weights[t]*generation[r,t]*co2_factors[r] for r in 1:G, t in 1:T)
    output["Load shedding"] = value.(ED[:total_nse])
    output["Net Flow"] = value.(ED[:flow])
    if crm_flag
        output["CRM dual"] = dual.(ED[:cCapacityResMargin]).*settings["Scaling factor cost"]./t_weights
        output["CRM slack"] = value.(ED[:vCRMSlack])
    end
    sd = value.(ED[:served_demand]).*settings["Scaling factor demand"]
    tb = value.(ED[:total_benefit]).*(settings["Scaling factor cost"]*settings["Scaling factor demand"])
    #output["Consumer surplus"] = sum(t_weights[t]*(tb[t] - output["Power price"][t]*sum(sd[seg,t] for seg in 1:nse_segs)) for t in 1:T)
    # Time-series outputs are large; only materialize when full scenario writes are requested.
    if settings["Write all scenarios flag"]
        output["Generation"] = generation
        if storage_flag
            output["Storage charging"] = value.(ED[:z_ch])
            output["Storage discharging"] = value.(ED[:z_dch])
        end
    end

    return output
end