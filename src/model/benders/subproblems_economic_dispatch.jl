export run_economic_dispatch, build_subproblems, set_capacity_parameters!, run_subproblem
using Distributed, DistributedArrays, JuMP, Gurobi
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
    if settings["Parallel flag"]
        @info "Building subproblems in parallel using $(Threads.nthreads()) workers..."
        Tot_scen = S*F*K
        SPs_flat = reshape(SPs, Tot_scen, 1)
        P = collect(1:Tot_scen)
        scenarios = [[div(p-1, F*K) + 1, div(mod(p-1, F*K), K) + 1, mod(p-1, K) + 1] for p in P]
        mapping = Dict{Int, Vector{Int}}(p => scenarios[p] for p in P)
        Threads.@threads for p in eachindex(SPs_flat)
            SPs_flat[p] = build_subproblem(inputs, settings, mapping[p]) 
        end
        SPs = reshape(SPs_flat, S, F, K)
    else
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

    demand_shift_weather_scenario = inputs["Demand shift weather"][:,k]
    variable_costs_scenario = inputs["Variable costs"][:,f]
    availability_scenario = inputs["Generation availability"][:,:,k]
    period_weights_scenario = inputs["Period weights"][:,k]
    demand_adder_scenario = inputs["Demand adders"][s]

     # Numerical settings
    scaling_factor_demand = settings["Scaling factor demand"]
    scaling_factor_cost = settings["Scaling factor cost"] 

    # Model
    ED = Model(Gurobi.Optimizer)
    ED.ext[:Demand_shift_weather_scenario] = demand_shift_weather_scenario
    ED.ext[:Variable_costs_scenario] = variable_costs_scenario
    ED.ext[:Availability_scenario] = availability_scenario
    ED.ext[:Period_weights_scenario] = period_weights_scenario
    ED.ext[:Demand_adder_scenario] = demand_adder_scenario
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
    @variable(ED, x[r in 1:R] in Parameter(0.0)) # Capacity, MW
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

    # ~~~
    # Objective function
    # ~~~ 

    @objective(ED, Min, sum(t_weights[t]*g[r,t]*cost_var[r] for r in 1:G, t in 1:T) + sum(t_weights[t]*cost_nse[t] for t in 1:T))

    return ED
end


function set_capacity_parameters!(SPs::Array{Model,3}, capacity::Vector{Float64})
    for SP in SPs
        set_parameter_value.(SP[:x], capacity)
    end
end

function run_subproblem(ED::Model, inputs, settings)
    output = Dict{String, Any}()

    optimize!(ED)
    if termination_status(ED) != MOI.OPTIMAL
        @warn("Subproblem did not solve to optimality. Status: $(termination_status(ED))")
        output["coeff"] = 0
    else
        output["coeff"] = 1
    end
    
    output_full = write_subproblem_results(ED,inputs, settings)
    output = merge(output, output_full)

    return output
end

function run_all_subproblems(SPs::Array{Model,3}, inputs, settings)

    sp_results = Array{Dict{String, Any},3}(undef, size(SPs)...)
    if settings["Parallel flag"]
        SPs_flat = reshape(SPs, :, 1)
        sp_results_flat = reshape(sp_results, :, 1)
        
        Threads.@threads for p in eachindex(SPs_flat)
            sp_results_flat[p] = run_subproblem(SPs_flat[p], inputs, settings)
        end
        sp_results = reshape(sp_results_flat, size(SPs))
    else
        for s in 1:size(SPs,1), f in 1:size(SPs,2), k in 1:size(SPs,3)
            sp_results[s,f,k] = run_subproblem(SPs[s,f,k], inputs, settings) 
        end
    end
    return sp_results
end

function write_subproblem_results(ED::Model, inputs, settings)
    output = Dict{String, Any}()

    
    availability_scenario = ED.ext[:Availability_scenario]
    period_weights_scenario = ED.ext[:Period_weights_scenario]

    # Settings
    storage_flag = inputs["Storage flag"]

    # Parameters
    co2_factors = inputs["CO2 emission intensities"]
    
    # Representative periods
    t_weights = period_weights_scenario


    # Sets
    T = inputs["Number of periods"]
    G = inputs["Number of generation resources"]
    O =  inputs["Number of storage resources"]
    storage_flag = O > 0


    max_nse = settings["Max nse"] # percentage of demand that can be curtailed in each segments
    nse_segs = length(max_nse)
    
    
    output["SP dual"] = dual.(ParameterRef.(ED[:x])).*settings["Scaling factor cost"]
    output["SP objective"] = objective_value(ED).*(settings["Scaling factor cost"]*settings["Scaling factor demand"])
    output["Power price"] = dual.(ED[:power_balance]).*settings["Scaling factor cost"]./t_weights
    nse_dual = dual.(ED[:max_load_shedding]).*settings["Scaling factor cost"]
    output["Max NSE dual"] = sum(max_nse[seg]*ED[:demand][t]*nse_dual[seg,t] for seg in 1:nse_segs, t in 1:T)
    generation = value.(ED[:g])
    output["Emissions"] = sum(t_weights[t]*generation[r,t]*co2_factors[r] for r in 1:G, t in 1:T)
    output["Load shedding"] = value.(ED[:total_nse])
    sd = value.(ED[:served_demand]).*settings["Scaling factor demand"]
    tb = value.(ED[:total_benefit]).*(settings["Scaling factor cost"]*settings["Scaling factor demand"])
    output["Consumer surplus"] = sum(t_weights[t]*(tb[t] - output["Power price"][t]*sum(sd[seg,t] for seg in 1:nse_segs)) for t in 1:T)
    output["Generation"] = generation
    if storage_flag
        output["Storage charging"] = value.(ED[:z_ch])
        output["Storage discharging"] = value.(ED[:z_dch])
    end

    return output
end