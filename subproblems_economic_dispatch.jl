using JuMP
using CSV
using DataFrames
using Gurobi

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
    @variable(ED, x[r in 1:R] >= 0)
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
    end

    # Load shedding
    # Maximum non served energy (primal feasibility)
    @constraint(ED, max_load_shedding[seg in 1:nse_segs, t in 1:T, z in 1:Z], max_nse[seg]*demand[t,z] - y[seg,t,z] >= 0)
    @expression(ED, total_nse[t in 1:T, z in 1:Z], sum(y[seg,t,z] for seg in 1:nse_segs))
    # Consumer surplus including demand shock for reporting purposes
    @expression(ED, served_demand[seg in 1:nse_segs, t in 1:T, z in 1:Z], max_nse[seg]*demand[t,z] - y[seg,t,z])
    @expression(ED, total_benefit[t in 1:T, z in 1:Z], price_nse[1]*demand_base[t,z] + sum(price_nse[seg]*served_demand[seg,t,z] for seg in 2:nse_segs))

    # Capacity limit on generation
    @constraint(ED, capacity_limit[r in 1:G, t in 1:T], g[r,t] <= x[r]*availability[t,r])
    if Z > 1
        @variable(ED, flow[t in 1:T, l in 1:L])
        @constraint(ED, transmission_limit[t in 1:T, l in 1:L], flow[t,l] <= inputs["Line capacities"][l]/scaling_factor_demand)
        @constraint(ED, transmission_limit_negative[t in 1:T, l in 1:L], flow[t,l] >= -inputs["Line capacities"][l]/scaling_factor_demand)
    end

    # Power balance constraint (per zone)
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

    @constraint(ED, planning_problem_capacity[r in 1:R], x[r] == capacity[r])

    # ~~~
    # Objective function
    # ~~~ 

    @objective(ED, Min, sum(t_weights[t]*g[r,t]*cost_var[r] for r in 1:G, t in 1:T) + sum(t_weights[t]*cost_nse[t,z] for t in 1:T, z in 1:Z))

    optimize!(ED)

    # Write outputs
    output = Dict{String, Any}()
    output["SP dual"] = dual.(planning_problem_capacity).*scaling_factor_cost
    output["SP objective"] = objective_value(ED).*(scaling_factor_cost*scaling_factor_demand)
    output["Power price"] = dual.(power_balance).*scaling_factor_cost./t_weights
    nse_dual = dual.(max_load_shedding).*scaling_factor_cost
    generation = value.(g)
    output["Emissions"] = sum(t_weights[t]*generation[r,t]*co2_factors[r] for r in 1:G, t in 1:T)
    output["Load shedding"] = value.(total_nse)

    return output
end
