export run_economic_dispatch, build_subproblems, set_capacity_parameters!, run_subproblem
using Distributed, DistributedArrays, JuMP, Gurobi

const WORKER_SP_CACHE = Dict{Tuple{Int, Int, Int}, Model}()
const WORKER_DATA_CACHE = Dict{Symbol, Any}()
const WORKER_MGCA_CACHE = Dict{Symbol, Any}()
const WORKER_MODEL_CACHE = Dict{Symbol, Any}()
const MASTER_WORKER_SP_KEYS = Dict{Int, Set{Tuple{Int, Int, Int}}}()
const MASTER_SCENARIO_TO_WORKER = Dict{Tuple{Int, Int, Int}, Int}()
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

function _worker_log_memory!(label::String, scenario_key::Tuple{Int, Int, Int}; payload_size_mb::Union{Nothing, Float64}=nothing)
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

function _build_and_cache_subproblem!(scenario_key::Tuple{Int, Int, Int})
    SP = build_subproblem(WORKER_DATA_CACHE[:inputs], WORKER_DATA_CACHE[:settings], collect(scenario_key))
    WORKER_SP_CACHE[scenario_key] = SP
    return nothing
end

function _run_cached_subproblem!(scenario_key::Tuple{Int, Int, Int}, capacity::Vector{Float64}, capacity_line::Vector{Float64}, minimal_payload::Bool)
    SP = WORKER_SP_CACHE[scenario_key]
    set_parameter_value.(SP[:x], capacity)
    set_parameter_value.(SP[:x_line], capacity_line)
    result = run_subproblem(SP, WORKER_DATA_CACHE[:inputs], WORKER_DATA_CACHE[:settings]; minimal_payload=minimal_payload)
    return result
end

function _clear_worker_mgca_cache!()
    empty!(WORKER_MGCA_CACHE)
    GC.gc(false)
    return nothing
end

function _set_worker_mgca_problem_data!(points, solver)
    WORKER_MGCA_CACHE[:points] = points
    WORKER_MGCA_CACHE[:solver] = solver

    return nothing
end

function _build_and_cache_mgca_problem!()
    points = WORKER_MGCA_CACHE[:points]
    solver = WORKER_MGCA_CACHE[:solver]

    rows, cols = size(points)
    bounded_points = max.(points, 0.0)
    model = Model()
    if solver == "HiGHS"
        set_optimizer(model, Ipopt.Optimizer)
    elseif solver == "Gurobi"
        set_optimizer(model, Gurobi.Optimizer)
    end
    @variable(model, l[1:rows] >= 0)
    @constraint(model, sum(l[i] for i in 1:rows) == 1)
    @variable(model, x[1:cols])
    @constraint(model, x == bounded_points' * l)
    @constraint(model, l .<= 0.95)
    @objective(model, Min, 0)
    set_silent(model)

    WORKER_MGCA_CACHE[:model] = model
    return nothing
end

function _run_cached_mgca_sample!(target::Vector{Float64})
    model = WORKER_MGCA_CACHE[:model]
    @objective(model, Min, sum((model[:x][i] - target[i])^2 for i in 1:length(model[:x])))
    optimize!(model)
    return max.(value.(model[:x]), 0.0)
end

##### TODO: Make better sampling methods - specifically want to try centroidal voronoi tessellations (CVT) which are designed to produce well-distributed samples in high-dimensional spaces. The current method of sampling random points and projecting them onto the convex hull of the original points is a simple approach but may not yield the best distribution of samples, especially as the dimensionality increases. CVT can help in generating samples that are more representative of the underlying distribution of the data, which could lead to better performance in downstream tasks that rely on these samples.
# try on both weights and original points
function cvt_sampler(points, num_samples; max_iters = 100)
    k = num_samples
    (n_points, n_vars) = size(points)
    num_initials = k * 10 
    max_point = maximum(points, dims=1)[:]
    initial_centers = Matrix([rand(length(max_point)) .* max_point for _ in 1:num_initials])
    R = kmeans(initial_centers', k; maxiter=max_iters)
    samples = collect(col for col in eachcol(R.centers))
    return samples
end

# One hit-and-run step in weight-space on the polytope l_i in [0, ub], sum(l) = 1 (the
# same simplex-with-cap the MGCA projection uses to express the convex hull of `points`).
# Directions are drawn from the sum-zero hyperplane so every step stays on sum(l) = 1;
# for rows == 2 this polytope is just the segment between (ub, 1-ub) and (1-ub, ub), so
# the walk correctly degenerates to a 1-D bounce rather than exploring a higher-dim body.


# Hit and run as parameterizedr right now is very bad for our purposes.
function _hit_and_run_step(l::Vector{Float64}, rng::AbstractRNG; ub::Float64 = 0.95, max_direction_retries::Int = 50)
    rows = length(l)
    for _ in 1:max_direction_retries
        d = randn(rng, rows)
        d .-= sum(d) / rows
        dn = norm(d)
        dn < 1e-10 && continue
        d ./= dn

        t_min, t_max = -Inf, Inf
        for i in 1:rows
            if d[i] > 1e-12
                t_min = max(t_min, -l[i] / d[i])
                t_max = min(t_max, (ub - l[i]) / d[i])
            elseif d[i] < -1e-12
                t_min = max(t_min, (ub - l[i]) / d[i])
                t_max = min(t_max, -l[i] / d[i])
            end
        end

        if t_max - t_min < 1e-9
            continue
        end

        t = t_min + rand(rng) * (t_max - t_min)
        l_new = clamp.(l .+ t .* d, 0.0, ub)
        if abs(sum(l_new) - 1.0) > 1e-8
            l_new ./= sum(l_new)
        end
        return l_new
    end
    return l
end

# Runs a single hit-and-run chain of length rows in l-space, discarding burn_in steps and
# keeping every thin-th step thereafter. Returns n_samples weight vectors summing to 1.
function _run_hit_and_run_chain(rows::Int, n_samples::Int; rng::AbstractRNG, burn_in::Int, thin::Int, ub::Float64 = 0.95)
    rows == 1 && error("Hit-and-run sampling requires rows > 1 (got rows=1): the simplex-with-cap polytope l_1=1, l_1<=$(ub) is empty.")

    l = fill(1.0 / rows, rows)
    for _ in 1:burn_in
        l = _hit_and_run_step(l, rng; ub = ub)
    end

    samples = Vector{Vector{Float64}}(undef, n_samples)
    for idx in 1:n_samples
        for _ in 1:thin
            l = _hit_and_run_step(l, rng; ub = ub)
        end
        samples[idx] = l
    end
    return samples
end

# Hit-and-run MCMC sampler over the convex hull of `points`, walked in weight-space (see
# _hit_and_run_step) then mapped to capacity-space via x = bounded_points' * l. Unlike the
# CVT/Random methods this needs no solver at all - each returned sample is a pure convex
# combination of `points`. `rng` defaults to a seed derived from myid(), so dispatching
# this to different workers via remotecall/@fetchfrom gives each worker an independent,
# reproducible chain for free.
function hit_and_run_sampler(points, num_samples, settings; rng::Union{Nothing, AbstractRNG} = nothing)
    rows, cols = size(points)
    bounded_points = max.(points, 0.0)

    burn_in = get(settings, "Hit and Run Burn-in", max(500, 20 * rows))
    thin = get(settings, "Hit and Run Thin", max(10, 2 * rows))
    chain_rng = something(rng, MersenneTwister(get(settings, "Seed", 1234) + myid()))

    l_samples = _run_hit_and_run_chain(rows, num_samples; rng = chain_rng, burn_in = burn_in, thin = thin)
    return [bounded_points' * l for l in l_samples]
end

function sample_interior_distributed(points, num_samples, settings)
    samples = Vector{Vector{Float64}}(undef, num_samples)
    solver = settings["Solver"]
    sample_method = settings["Sample Method"]
    if nprocs() > 1
        pids = workers()
        if sample_method == "HitAndRun"
            @info "Sampling using hit-and-run MCMC method with $num_samples samples across $(length(pids)) workers..."
            n_workers = length(pids)
            per_worker = ceil(Int, num_samples / n_workers)
            worker_batches = Vector{Vector{Vector{Float64}}}(undef, n_workers)
            @sync for (w, pid) in enumerate(pids)
                @async worker_batches[w] = @fetchfrom pid hit_and_run_sampler(points, per_worker, settings)
            end
            samples = vcat(worker_batches...)[1:num_samples]
        elseif sample_method in ("CVT", "Random")
            if sample_method == "CVT"
                @info "Sampling using CVT method with $num_samples samples across $(length(pids)) workers..."
                sample_targets = cvt_sampler(points, num_samples)
            else
                @info "Sampling using random projection method with $num_samples samples across $(length(pids)) workers..."
                max_point = maximum(points, dims=1)[:]
                sample_targets = [rand(length(max_point)) .* max_point for _ in 1:num_samples]
            end

            @sync for pid in pids
                @async begin
                    remotecall_wait(_clear_worker_mgca_cache!, pid)
                    remotecall_wait(_set_worker_mgca_problem_data!, pid, points, solver)
                    remotecall_wait(_build_and_cache_mgca_problem!, pid)
                end
            end

            @sync for (idx, target) in enumerate(sample_targets)
                @async begin
                    pid = pids[mod1(idx, length(pids))]
                    samples[idx] = @fetchfrom pid _run_cached_mgca_sample!(target)
                end
            end
        else
            error("Unsupported sample method: $sample_method. Supported methods are: CVT, Random, HitAndRun")
        end
        @everywhere GC.gc()
    else
        if sample_method == "HitAndRun"
            @info "Sampling using hit-and-run MCMC method with $num_samples samples (serial, single chain)..."
            samples = hit_and_run_sampler(points, num_samples, settings)
        elseif sample_method in ("CVT", "Random")
            model = Model()
            if solver == "HiGHS"
                set_optimizer(model, Ipopt.Optimizer)
            elseif solver == "Gurobi"
                set_optimizer(model, Gurobi.Optimizer)
            end
            rows, cols = size(points)
            bounded_points = max.(points, 0.0)
            @variable(model, l[1:rows] >= 0)
            @constraint(model, sum(l[i] for i in 1:rows) == 1)
            @variable(model, x[1:cols])
            @constraint(model, x == bounded_points' * l)
            @constraint(model, l .<= 0.95)
            @objective(model, Min, 0)
            set_silent(model)

            if sample_method == "CVT"
                @info "Sampling using CVT method with $num_samples samples (serial)..."
                sample_targets = cvt_sampler(points, num_samples)
            else
                @info "Sampling using random projection method with $num_samples samples (serial)..."
                max_point = maximum(points, dims=1)[:]
                sample_targets = [rand(length(max_point)) .* max_point for _ in 1:num_samples]
            end

            for (idx, target) in enumerate(sample_targets)
                @objective(model, Min, sum((model[:x][i] - target[i])^2 for i in 1:length(model[:x])))
                optimize!(model)
                samples[idx] = max.(value.(model[:x]), 0.0)
            end
        else
            error("Unsupported sample method: $sample_method. Supported methods are: CVT, Random, HitAndRun")
        end
    end

    return samples
end

function run_economic_dispatch(inputs, settings, capacity, demand_shift_weather_scenario, variable_costs_scenario, availability_scenario, period_weights_scenario, demand_adder_scenario)

    # Check for negative capacity and set to 0
    for r in 1:length(capacity)
        if capacity[r] < 0
            capacity[r] = 0.0
        end
    end

    # Numerical settings
    scaling_factor_demand = settings["Scaling factor demand"]
    scaling_factor_cost = settings["Scaling factor cost"] 

    # Model
    ED = Model(Gurobi.Optimizer)
    set_optimizer_attribute(ED, "OutputFlag", 0)
    set_optimizer_attribute(ED, "Crossover", 1)
    set_optimizer_attribute(ED, "Method", 2)
    set_optimizer_attribute(ED, "BarConvTol", 1e-5)
    set_optimizer_attribute(ED, "OptimalityTol", 1e-5)
    set_optimizer_attribute(ED, "FeasibilityTol", 1e-5)
    set_silent(ED)

    # ~~~~
    # Load inputs
    # ~~~

    # Settings
    storage_flag = inputs["Storage flag"]

    # Parameters
    availability = availability_scenario
    co2_factors = inputs["CO2 emission intensities"]
    
    # Storage specific
    if storage_flag
        N = inputs["Storage power to energy ratio"]
        # Efficiency
        F_dch = inputs["Storage discharge loss"]
        F_ch = inputs["Storage charge loss"]

        first_periods = inputs["Set of first periods"]
        last_periods = inputs["Set of last periods"]
    end

    # Representative periods
    n_periods = inputs["Number of representative periods"]
    t_weights = period_weights_scenario

    # ~~~
    # Build model
    # ~~~

    # Sets
    T = inputs["Number of periods"]
    G = inputs["Number of generation resources"]
    O =  inputs["Number of storage resources"]
    storage_flag = O > 0
    R = G + O

    # ~~~
    # Parameters
    # ~~~
    
    # Scale parameters
    
    demand_base = inputs["Demand"]./scaling_factor_demand
    demand_add = demand_adder_scenario./scaling_factor_demand
    demand_shift_weather = demand_shift_weather_scenario./scaling_factor_demand
    cost_var = variable_costs_scenario./scaling_factor_cost
    
    # Demand
    @expression(ED, demand[t in 1:T], demand_base[t] + demand_shift_weather[t] + demand_add)
    max_nse = settings["Max nse"] # percentage of demand that can be curtailed in each segments
    nse_segs = length(max_nse)
    price_nse = settings["Cost of nonserved energy"]./scaling_factor_cost
    
    # ~~~
    # Model formulation
    # ~~~

    # Generation
    @variable(ED, g[r in 1:G, t in 1:T] >= 0) # MWh
    @variable(ED, x[r in 1:R] >= 0)
    # Non-servED energy
    @variable(ED, y[seg in 1:nse_segs, t in 1:T] >= 0) # $/MWh
    @expression(ED, cost_nse[t in 1:T], sum(price_nse[seg]*y[seg,t] for seg in 1:nse_segs))

    if storage_flag
        @variable(ED, z_dch[t in 1:T] >= 0)
        @variable(ED, z_ch[t in 1:T] >= 0)
        # State of charge
        @variable(ED, e[t in 1:T] >= 0) 
        # Primal feasibility
        
        for i in 1:n_periods
            # Wrap first and last periods - start
            @constraints(ED, 
            begin 
                e[first_periods[i]] == e[last_periods[i]] - (1/F_dch)*z_dch[first_periods[i]] + F_ch*z_ch[first_periods[i]]
            end)
            @constraints(ED, 
            begin z_dch[first_periods[i]] <= e[last_periods[i]]
            end)
            # Wrap first and last periods - end
            @constraints(ED, 
            begin [t in first_periods[i]+1:last_periods[i]], e[t] == e[t-1] - (1/F_dch)*z_dch[t] + F_ch*z_ch[t] end)
            @constraints(ED, begin [t in first_periods[i]+1:last_periods[i]], z_dch[t] <= e[t-1] end)
        end

        
        # Energy balance for the remaining periods
        @constraint(ED, energy_limit[t in 1:T], e[t] <= (1/N)*x[G+1])
        @constraint(ED, charge_limit_rate[t in 1:T], z_ch[t] <= x[G+1])
        
        @constraint(ED, discharge_limit_rate[t in 1:T], z_dch[t] <= x[G+1])
        @constraint(ED, charge_discharge_balance[t in 1:T], z_dch[t] + z_ch[t] <= x[G+1])
    end

    # Load shedding
    # Maximum non served energy (primal feasibility)
    @constraint(ED, max_load_shedding[seg in 1:nse_segs,t in 1:T], max_nse[seg]*demand[t] - y[seg,t] >= 0)
    @expression(ED, total_nse[t in 1:T], sum(y[seg,t] for seg in 1:nse_segs))
    # Consumer surplus including demand shock for reporting purposes
    @expression(ED, served_demand[seg in 1:nse_segs, t in 1:T], max_nse[seg]*demand[t] - y[seg,t])
    @expression(ED, total_benefit[t in 1:T], price_nse[1]*demand_base[t] + sum(price_nse[seg]*served_demand[seg,t] for seg in 2:nse_segs))

    # Capacity limit on generation
    @constraint(ED, capacity_limit[r in 1:G, t in 1:T], g[r,t] <= x[r]*availability[t,r])

    # Power balance constraint
    if storage_flag
        @constraint(ED, power_balance[t in 1:T], sum(g[r,t] for r in 1:G) + total_nse[t] + z_dch[t] - z_ch[t] - demand[t] == 0)
    else
        @constraint(ED, power_balance[t in 1:T], sum(g[r,t] for r in 1:G) + total_nse[t] - demand[t] == 0)
    end

    @constraint(ED, planning_problem_capacity[r in 1:R], x[r] == capacity[r])

    # ~~~
    # Objective function
    # ~~~ 

    @objective(ED, Min, sum(t_weights[t]*g[r,t]*cost_var[r] for r in 1:G, t in 1:T) + sum(t_weights[t]*cost_nse[t] for t in 1:T))

    optimize!(ED)

    # Write outputs
    output = Dict{String, Any}()
    output["SP dual"] = dual.(planning_problem_capacity).*scaling_factor_cost
    output["SP objective"] = objective_value(ED).*(scaling_factor_cost*scaling_factor_demand)
    output["Power price"] = dual.(power_balance).*scaling_factor_cost./t_weights
    nse_dual = dual.(max_load_shedding).*scaling_factor_cost
    output["Max NSE dual"] = sum(max_nse[seg]*demand[t]*nse_dual[seg,t] for seg in 1:nse_segs, t in 1:T)
    generation = value.(g)
    output["Emissions"] = sum(t_weights[t]*generation[r,t]*co2_factors[r] for r in 1:G, t in 1:T)
    output["Load shedding"] = value.(total_nse)
    sd = value.(served_demand).*scaling_factor_demand
    tb = value.(total_benefit).*(scaling_factor_cost*scaling_factor_demand)
    output["Consumer surplus"] = sum(t_weights[t]*(tb[t] - output["Power price"][t]*sum(sd[seg,t] for seg in 1:nse_segs)) for t in 1:T)

    return output
end

function build_all_subproblems(inputs, settings)
    S = inputs["Number of demand scenarios"]
    F = inputs["Number of gas price scenarios"]
    K = inputs["Number of weather scenarios"]
    SPs = Array{Model,3}(undef, S, F, K)
    if settings["Parallel flag"] && nworkers() > 0
        @info "Building subproblems in parallel using $(nworkers()) distributed workers..."
        pids = workers()
        empty!(MASTER_SCENARIO_TO_WORKER)
        @sync for pid in pids
            @async begin
                remotecall_wait(_clear_worker_subproblem_cache!, pid)
                remotecall_wait(_set_worker_problem_data!, pid, inputs, settings)
                worker_set = get!(MASTER_WORKER_SP_KEYS, pid, Set{Tuple{Int, Int, Int}}())
                empty!(worker_set)
            end
        end
        scenario_keys = [(s, f, k) for s in 1:S, f in 1:F, k in 1:K]
        @sync for (idx, scenario_key) in enumerate(scenario_keys)
            @async begin
                s, f, k = scenario_key
                SPs[s,f,k] = build_subproblem(inputs, settings, [s, f, k])
                pid = pids[mod1(idx, length(pids))]
                remotecall_wait(_build_and_cache_subproblem!, pid, scenario_key)
                push!(MASTER_WORKER_SP_KEYS[pid], scenario_key)
                MASTER_SCENARIO_TO_WORKER[scenario_key] = pid
            end
        end
    else
        if settings["Parallel flag"] && nworkers() == 0
            @warn("Parallel flag is true but no distributed workers are available. Falling back to serial subproblem build.")
        end
         for s in 1:S, f in 1:F, k in 1:K
            SPs[s,f,k] = build_subproblem(inputs, settings, [s,f,k])
        end
    end
    return SPs
end

function build_subproblem(inputs, settings, scenario_index::AbstractVector)
    # Set Scenario index
    s = scenario_index[1]
    f = scenario_index[2]
    k = scenario_index[3]

    demand_shift_weather_scenario = inputs["Demand shift weather"][:,:,k]
    variable_costs_scenario = inputs["Variable costs"][:,f]
    availability_scenario = inputs["Generation availability"][:,:,k]
    period_weights_scenario = inputs["Period weights"][:,k]
    demand_adder_scenario = inputs["Demand adders"][s]

     # Numerical settings
    scaling_factor_demand = settings["Scaling factor demand"]
    scaling_factor_cost = settings["Scaling factor cost"] 

    # Model
    ED = Model()
    ED.ext[:Demand_shift_weather_scenario] = demand_shift_weather_scenario
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

        first_periods = inputs["Set of first periods"]
        last_periods = inputs["Set of last periods"]
    end

    # Representative periods
    n_periods = inputs["Number of representative periods"]
    t_weights = period_weights_scenario
    zone_map = inputs["Zone map"]
    res_zone_map = inputs["Resource zone map"]

    # ~~~
    # Build model
    # ~~~

    # Sets
    T = inputs["Number of periods"]
    G = inputs["Number of generation resources"]
    O =  inputs["Number of storage resources"]
    storage_flag = O > 0
    R = G + O
    L = inputs["Number of lines"]
    Z = inputs["Number of zones"]

    # ~~~
    # Parameters
    # ~~~
    
    # Scale parameters
    
    demand_base = inputs["Demand"]./scaling_factor_demand
    demand_add = demand_adder_scenario./scaling_factor_demand
    demand_shift_weather = demand_shift_weather_scenario./scaling_factor_demand
    cost_var = variable_costs_scenario./scaling_factor_cost
    # Demand
    @expression(ED, demand[t in 1:T, z in 1:Z], demand_base[t, z] + demand_shift_weather[t, z] + demand_add)
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
        @variable(ED, z_dch[t in 1:T] >= 0)
        @variable(ED, z_ch[t in 1:T] >= 0)
        # State of charge
        @variable(ED, e[t in 1:T] >= 0) 
        # Primal feasibility
        
        for i in 1:n_periods
            # Wrap first and last periods - start
            @constraints(ED, 
            begin 
                e[first_periods[i]] == e[last_periods[i]] - (1/F_dch)*z_dch[first_periods[i]] + F_ch*z_ch[first_periods[i]]
            end)
            @constraints(ED, 
            begin z_dch[first_periods[i]] <= e[last_periods[i]]
            end)
            # Wrap first and last periods - end
            @constraints(ED, 
            begin [t in first_periods[i]+1:last_periods[i]], e[t] == e[t-1] - (1/F_dch)*z_dch[t] + F_ch*z_ch[t] end)
            @constraints(ED, begin [t in first_periods[i]+1:last_periods[i]], z_dch[t] <= e[t-1] end)
        end

        
        # Energy balance for the remaining periods
        @constraint(ED, energy_limit[t in 1:T], e[t] <= (1/N)*x[G+1])
        @constraint(ED, charge_limit_rate[t in 1:T], z_ch[t] <= x[G+1])
        
        @constraint(ED, discharge_limit_rate[t in 1:T], z_dch[t] <= x[G+1])
        @constraint(ED, charge_discharge_balance[t in 1:T], z_dch[t] + z_ch[t] <= x[G+1])
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
    if Z > 1
        @variable(ED, flow[t in 1:T, l in 1:L])
        @constraint(ED, transmission_limit[t in 1:T, l in 1:L], flow[t,l] <= x_line[l]/scaling_factor_demand)
        @constraint(ED, transmission_limit_negative[t in 1:T, l in 1:L], flow[t,l] >= -x_line[l]/scaling_factor_demand)
    end

    # Power balance constraint
    @expression(ED, supply[t in 1:T, z in 1:Z], sum(g[r,t]*res_zone_map[r,z] for r in 1:G))
    if storage_flag
        for t in 1:T, z in 1:Z
            add_to_expression!(supply[t,z], res_zone_map[G+1,z]*(z_dch[t] - z_ch[t]))
        end
    end
    if Z > 1
        for t in 1:T, z in 1:Z
            add_to_expression!(supply[t,z], sum(flow[t,l]*zone_map[l,z] for l in 1:L))
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
            sum(derate_gen[r]*x[r]*availability[t,r]*res_zone_map[r,z] for r in 1:G)
            + sum(y[seg,t,z] for seg in 2:nse_segs))
        if storage_flag
            for t in 1:T, z in 1:Z
                add_to_expression!(eCapResMargin[t,z], res_zone_map[G+1,z]*(z_dch[t] - z_ch[t]))
            end
        end
        if Z > 1
            for t in 1:T, z in 1:Z
                add_to_expression!(eCapResMargin[t,z], sum(flow[t,l]*zone_map[l,z] for l in 1:L))
            end
        end

        @variable(ED, vCRMSlack[t in 1:T, z in 1:Z] >= 0)
        @constraint(ED, cCapacityResMargin[t in 1:T, z in 1:Z],
            eCapResMargin[t,z] + vCRMSlack[t,z] >= (1 + crm_margin[z])*demand[t,z])
        @expression(ED, crm_slack_cost, sum(t_weights[t]*crm_price_cap[z]*vCRMSlack[t,z] for t in 1:T, z in 1:Z))
    else
        @expression(ED, crm_slack_cost, AffExpr(0.0))
    end


    @expression(ED, cost_by_zone[z in 1:Z], sum(g[r,t]*res_zone_map[r,z]*cost_var[r] for r in 1:G, t in 1:T) + sum(cost_nse[t,z] for t in 1:T))
    @expression(ED, emissions_by_zone[z in 1:Z], sum(t_weights[t]*sum(g[r,t]*res_zone_map[r,z]*co2_factors[r] for r in 1:G) for t in 1:T))
    # ~~~
    # Objective function
    # ~~~

    @objective(ED, Min, sum(t_weights[t]*g[r,t]*cost_var[r] for r in 1:G, t in 1:T) + sum(t_weights[t]*cost_nse[t,z] for t in 1:T, z in 1:Z) + crm_slack_cost)

    return ED
end


function set_capacity_parameters!(SPs::Array{Model,3}, capacity::Vector{Float64}, capacity_line::Vector{Float64})
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

function run_all_subproblems(SPs::Array{Model,3}, inputs, settings, capacity::Vector{Float64}, capacity_line::Vector{Float64}; minimal_payload::Bool=get(settings, "Minimal worker payload flag", true))

    sp_results = Array{Dict{String, Any},3}(undef, size(SPs)...)
    if settings["Parallel flag"] && nworkers() > 0
        pids = workers()
        scenario_keys = [(s, f, k) for s in 1:size(SPs,1), f in 1:size(SPs,2), k in 1:size(SPs,3)]

        @sync for (idx, scenario_key) in enumerate(scenario_keys)
            @async begin
                if !haskey(MASTER_SCENARIO_TO_WORKER, scenario_key)
                    error("Subproblem $(scenario_key) has no worker assignment. Call build_all_subproblems() before run_all_subproblems().")
                end
                pid = MASTER_SCENARIO_TO_WORKER[scenario_key]
                if !(pid in pids)
                    error("Assigned worker $(pid) for subproblem $(scenario_key) is not active. Rebuild subproblems with build_all_subproblems().")
                end
                worker_keys = get(MASTER_WORKER_SP_KEYS, pid, Set{Tuple{Int, Int, Int}}())
                if !(scenario_key in worker_keys)
                    error("Subproblem $(scenario_key) is not cached on worker $(pid). Call build_all_subproblems() before run_all_subproblems().")
                end
                sp_results[scenario_key...] = @fetchfrom pid _run_cached_subproblem!(scenario_key, capacity, capacity_line, minimal_payload)
                
            end
        end
        @everywhere GC.gc()
    else
        if settings["Parallel flag"] && nworkers() == 0
            @warn("Parallel flag is true but no distributed workers are available. Falling back to serial subproblem solves.")
        end
        set_capacity_parameters!(SPs, capacity, capacity_line)
        for s in 1:size(SPs,1), f in 1:size(SPs,2), k in 1:size(SPs,3)
            sp_results[s,f,k] = run_subproblem(SPs[s,f,k], inputs, settings; minimal_payload=minimal_payload) 
        end
    end
    return sp_results
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


    # Sets
    T = inputs["Number of periods"]
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