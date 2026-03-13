using JuMP
using CSV
using DataFrames
using Gurobi
using Random
using LinearAlgebra

function optimization_model(inputs, settings)

    # Numerical settings
    scaling_factor_demand = 1 # to GWs
    scaling_factor_obj = 1 
    scaling_factor_cost = 1 # to $1,000,000s

    # Model
    gep = Model(Gurobi.Optimizer)

    set_optimizer_attribute(gep, "Method", 2)
    set_optimizer_attribute(gep, "Crossover", 0)
    set_optimizer_attribute(gep, "BarConvTol", 1e-2)
    
    # ~~~~
    # Load inputs
    # ~~~

    risk_aversion_flag = settings["Risk aversion flag"]
    risk_aversion_weight = settings["Risk aversion weight"]

    # Settings
    storage_flag = inputs["Storage flag"]
    itc_flag = settings["Investment tax credits flag"]
    ptc_flag = settings["Production tax credits flag"]
    co2_tax_flag = settings["Carbon tax flag"]
    cap_flag = settings["Carbon cap flag"]
    standard_flag = settings["Standard flag"]

    # General
    cost_inv = inputs["Investment costs"]
    cost_var = inputs["Variable costs"]
    co2_factors = inputs["CO2 emission intensities"]
    demand_orig_scale = inputs["Demand"]
    # Variability (1s for dispatchable and fraction for VRE)
    availability = inputs["Generation availability"]

    # Storage specific
    if storage_flag
        N = inputs["Storage power to energy ratio"]
        # Efficiency
        F_dch = inputs["Storage discharge loss"]
        F_ch = inputs["Storage charge loss"]
    end

    # Stochasticities
    P = inputs["Demand scenario probabilities"]
    P_f = inputs["Fuel price scenario probabilities"]

    # Representative periods
    n_periods = inputs["Number of representative periods"]
    t_weights = inputs["Period weights"]
    first_periods = inputs["Set of first periods"]
    last_periods = inputs["Set of last periods"]

    # ~~~
    # Build model
    # ~~~

    # Sets
    T = inputs["Number of periods"]
    S = size(P)[1] # number of demand scenarios
    F = size(P_f)[1] # number of fuel cost scenarios
    G = inputs["Number of generation resources"]
    O =  inputs["Number of storage resources"]
    storage_flag = O > 0
    R = G + O


    # ~~~
    # Parameters
    # ~~~
    # Policies
    
    if co2_tax_flag 
        co2_tax = settings["Carbon tax"]
    else
        co2_tax = 0
    end
    if itc_flag
        itc = settings["Tax credits"]
    else
        itc = zeros(R)
    end
    if ptc_flag
        ptc = settings["Production credits"]
    else
        ptc = zeros(R)
    end
    if standard_flag
        standard_perc = settings["Standard mandate"]
    else
        standard_perc = 0
    end

    if risk_aversion_flag
        # CVaR parameters
        Ψ = settings["Value-at-Risk percent"] 
        Ω = risk_aversion_weight
    else
        Ω = 1
    end

    # Scale parameters
    co2_tax = co2_tax/scaling_factor_cost
    demand = demand_orig_scale./scaling_factor_demand
    cost_inv = cost_inv./scaling_factor_cost
    cost_var = cost_var./scaling_factor_cost
    g_max = 50e3/scaling_factor_demand
    # Demand flexibility
    price_nse = [2e3,1800,1100,400]./scaling_factor_cost # originally $/MWh    
    max_nse = [1,0.04,0.024,0.003] # percentage of demand that can be curtailed in each segments
    NSE_SEGS = 4

    # ~~~
    # Model formulation
    # ~~~

    # Generation
    @variable(gep, 0 <= g[r in 1:G, t in 1:T, s in 1:S, f in 1:F] <= g_max) # MWh
    # Capacity
    @variable(gep, 0 <= x[r in 1:R]) # Capacity, MW
    # Non-served energy
    @variable(gep, y[seg in 1:NSE_SEGS, t in 1:T, s in 1:S, f in 1:F] >= 0) # $/MWh
    @expression(gep, eNSE_Cost[t in 1:T, s in 1:S, f in 1:F], sum(price_nse[seg]*y[seg,t,s,f] for seg in 1:NSE_SEGS))

    if risk_aversion_flag
        # Auxiliary varliables for CVaR
        @variable(gep, u[s in 1:S, f in 1:F] >= 0) # Gain/loss relative to VaR, $
        @variable(gep, ζ) # VaR variable, $
        @constraint(gep, cvar_tail[s in 1:S, f in 1:F], u[s,f] >= sum(t_weights[t]*g[r,t,s,f]*cost_var[r,f] for r in 1:G, t in 1:T) + sum(t_weights[t]*g[r,t,s,f]*co2_tax*co2_factors[r] for r in 1:G, t in 1:T) + sum(t_weights[t]*eNSE_Cost[t,s,f] for t in 1:T)  - ζ)
    end 

    if storage_flag
        @variable(gep, z_dch[t in 1:T, s in 1:S, f in 1:F] >= 0)
        @variable(gep, z_ch[t in 1:T, s in 1:S, f in 1:F] >= 0)
        # State of charge
        @variable(gep, e[t in 1:T, s in 1:S, f in 1:F] >= 0) 
        # Primal feasibility
        
        for i in 1:n_periods
            # Wrap first and last periods - start
            @constraints(gep, 
            begin 
                [s in 1:S, f in 1:F], e[first_periods[i],s,f] == e[last_periods[i],s,f] - (1/F_dch)*z_dch[first_periods[i],s,f] + F_ch*z_ch[first_periods[i],s,f]
            end)
            @constraints(gep, 
            begin [s in 1:S, f in 1:F], z_dch[first_periods[i],s,f] <= e[last_periods[i],s,f]
            end)
            # Wrap first and last periods - end
            @constraints(gep, 
            begin [t in first_periods[i]+1:last_periods[i], s in 1:S, f in 1:F], e[t,s,f] == e[t-1,s,f] - (1/F_dch)*z_dch[t,s,f] + F_ch*z_ch[t,s,f] end)
            @constraints(gep, begin [t in first_periods[i]+1:last_periods[i], s in 1:S, f in 1:F], z_dch[t,s,f] <= e[t-1,s,f] end)
        end

        
        # Energy balance for the remaining periods
        @constraint(gep, energy_limit[t in 1:T, s in 1:S, f in 1:F], e[t,s,f] <= (1/N)*x[G+1])
        @constraint(gep, charge_limit_rate[t in 1:T, s in 1:S, f in 1:F], z_ch[t,s,f] <= x[G+1])
        
        @constraint(gep, discharge_limit_rate[t in 1:T, s in 1:S, f in 1:F], z_dch[t,s,f] <= x[G+1])
        @constraint(gep, charge_discharge_balance[t in 1:T, s in 1:S, f in 1:F], z_dch[t,s,f] + z_ch[t,s,f] <= x[G+1])
    end

    # Load shedding
    # Maximum non served energy (primal feasibility)
    @constraint(gep, max_load_shedding[seg in 1:NSE_SEGS,t in 1:T, s in 1:S, f in 1:F], max_nse[seg]*demand[t,s] - y[seg,t,s,f] >= 0)
    @expression(gep,eNSE[t in 1:T,s in 1:S, f in 1:F],sum(y[seg,t,s,f] for seg in 1:NSE_SEGS))
    # Capacity limit on generation
    @expression(gep, capacity_limit_expression[r in 1:G, t in 1:T, s in 1:S, f in 1:F], x[r]*availability[t,r] - g[r,t,s,f])
    @constraint(gep, capacity_limit[r in 1:G, t in 1:T, s in 1:S, f in 1:F], capacity_limit_expression[r,t,s,f] >= 0)

    # Power balance constraint
    if storage_flag
        @constraint(gep, power_balance[t in 1:T, s in 1:S, f in 1:F], sum(g[r,t,s,f] for r in 1:G) + eNSE[t,s,f] + z_dch[t,s,f] - z_ch[t,s,f] - demand[t,s] == 0)
    else
        @constraint(gep, power_balance[t in 1:T, s in 1:S, f in 1:F], sum(g[r,t,s,f] for r in 1:G) + eNSE[t,s,f] - demand[t,s] == 0)
    end

    if standard_flag
        @constraint(gep, standard_constraint[s in 1:S, f in 1:F], sum(t_weights[t]*g[r,t,s,f] for r in 2:3, t in 1:T) >= standard_perc*sum(t_weights[t]*g[r,t,s,f] for r in 1:G, t in 1:T))
    end

    # Define system cost
    @expression(gep, op_cost[s in 1:S, f in 1:F], sum(t_weights[t]*g[r,t,s,f]*cost_var[r,f] for r in 1:G, t in 1:T) + sum(t_weights[t]*eNSE_Cost[t,s,f] for t in 1:T))
    @expression(gep, avg_op_cost[s in 1:S, f in 1:F], op_cost[s,f]/sum(t_weights[t]*demand[t,s] for t in 1:T))
    @expression(gep, inv_cost, sum(x[r]*cost_inv[r] for r in 1:R))
    @expression(gep, avg_inv_cost[s in 1:S], inv_cost/sum(t_weights[t]*demand[t,s] for t in 1:T))
    
    @expression(gep, system_cost, sum(x[r]*cost_inv[r] for r in 1:R) + sum(P[s]*P_f[f]*t_weights[t]*g[r,t,s,f]*cost_var[r,f] for r in 1:G, t in 1:T, s in 1:S, f in 1:F) + sum(P[s]*P_f[f]*t_weights[t]*eNSE_Cost[t,s,f] for t in 1:T, s in 1:S, f in 1:F))

    # Define emissions
    @expression(gep, emissions[s in 1:S, f in 1:F], sum(t_weights[t]*g[r,t,s,f]*co2_factors[r] for r in 1:G, t in 1:T))
    @expression(gep, emissions_exp, sum(P[s]*P_f[f]*emissions[s,f] for s in 1:S, f in 1:F))

    # Generation
    @expression(gep, eANN_GEN[r in 1:G, s in 1:S, f in 1:F], sum(g[r,t,s,f] for t in 1:T))

    
    # Cap constraint
    if cap_flag
        co2_cap = settings["Carbon cap"]
        @constraint(gep, emissions_cap[s in 1:S, f in 1:F], co2_cap >= emissions[s,f])
    else
        co2_cap = 0
    end

    # ~~~
    # Objective function
    # ~~~ 

    @expression(gep, sys_cost, sum(x[r]*cost_inv[r]*(1-itc[r]) for r in 1:R) + sum(P[s]*P_f[f]*t_weights[t]*g[r,t,s,f]*(cost_var[r,f]+co2_tax*co2_factors[r]) for r in 1:G, t in 1:T, s in 1:S, f in 1:F) + sum(P[s]*P_f[f]*t_weights[t]*eNSE_Cost[t,s,f] for t in 1:T, s in 1:S, f in 1:F))

    @expression(gep, risk_adjusted_sys_cost, sum(x[r]*cost_inv[r]*(1-itc[r]) for r in 1:R) + Ω*(sum(P[s]*P_f[f]*t_weights[t]*g[r,t,s,f]*(cost_var[r,f]+co2_tax*co2_factors[r]) for r in 1:G, t in 1:T, s in 1:S, f in 1:F) + sum(P[s]*P_f[f]*t_weights[t]*eNSE_Cost[t,s,f] for t in 1:T, s in 1:S, f in 1:F)) + (1-Ω)*(ζ + 1/Ψ*sum(P[s]*P_f[f]*u[s,f] for s in 1:S, f in 1:F)))

    if risk_aversion_flag
        @objective(gep, Min, risk_adjusted_sys_cost)
    else
        # Risk neutral
        @objective(gep, Min, sys_cost)
        
    end

    optimize!(gep)

    # Write outputs
    output = Dict{String, Any}()

    output["model"] = gep
        
    price = dual.(power_balance)
    price_shadow = shadow_price.(power_balance)  
    mu = dual.(capacity_limit)
    output["Capacity"] = value.(x)
    y = value.(eNSE)
    output["Load shedding"] = value.(eNSE)
    output["Generation"] = value.(g)
    output["Capacity dual"] = dual.(capacity_limit)
    output["Risk adjusted probabilities"] = zeros(S,F)
    output["Consumer surplus"] = zeros(S,F)
    if risk_aversion_flag
        theta = dual.(cvar_tail)/(1-Ω)
        for s in 1:S
            for f in 1:F
                output["Risk adjusted probabilities"][s,f] = theta[s,f]
            end
        end
    end


    if storage_flag
        output["Storage charging"] = value.(z_ch)
        output["Storage discharging"] = value.(z_dch)
    end
    output["Power balance dual"] = price
    output["Unweighted dual"] = scaling_factor_cost.*price./t_weights
    
    output["Objective function value"] = objective_value(gep)
    output["CVAR duals"] = dual.(cvar_tail)
    output["System cost"] = value.(system_cost)
    output["Operating cost"] = value.(op_cost)
    output["Operating cost average"] = value.(avg_op_cost)
    output["Investment cost"] = value.(inv_cost)
    output["Investment cost average"] = value.(avg_inv_cost)
    output["Scaling factor for costs"] = scaling_factor_cost
    output["Scaling factor for demand"] = scaling_factor_demand
    output["Emissions"] = value.(emissions)
    output["Emissions expected"] = value.(emissions_exp)
    # Record results from the linearization or not
    output["Piecewise linearized"] = false
    output["Risk sharing flag"] = true
    output["Optimization model flag"] = true
    if risk_aversion_flag
        output["CVaR Loss"] = value.(u)
        output["VaR"] = value.(ζ)
    else
        output["CVaR Loss"] = zeros(S,F)
        output["VaR"] = 0
    end

    if standard_flag
        output["Credit price"] = dual.(standard_constraint)*1e6
    else
        output["Credit price"] = 0
    end
    if cap_flag
        output["CO2 price"] = dual.(emissions_cap)*1e6
    else
        output["CO2 price"] = 0
    end
    
    return output
end


#### Build model ---- base stochastic model

function build_optimization_model(inputs, settings, objective_type::String, obj_weight::Float64)

    # Numerical settings
    scaling_factor_demand = 1 # to GWs
    scaling_factor_obj = 1 
    scaling_factor_cost = 1 # to $1,000,000s

    # Model
    gep = Model(Gurobi.Optimizer)

    set_optimizer_attribute(gep, "Method", 2)
    set_optimizer_attribute(gep, "Crossover", 0)
    set_optimizer_attribute(gep, "BarConvTol", 1e-5)
    set_optimizer_attribute(gep, "NumericFocus", 2)
    set_optimizer_attribute(gep, "FeasibilityTol", 1e-5)
    
    # ~~~~
    # Load inputs
    # ~~~

    risk_aversion_flag = objective_type
    risk_aversion_weight = obj_weight #(between 0 and 1)
    # Settings
    storage_flag = inputs["Storage flag"]
    itc_flag = settings["Investment tax credits flag"]
    ptc_flag = settings["Production tax credits flag"]
    co2_tax_flag = settings["Carbon tax flag"]
    cap_flag = settings["Carbon cap flag"]
    standard_flag = settings["Standard flag"]

    # General
    cost_inv = inputs["Investment costs"]
    cost_var = inputs["Variable costs"]
    co2_factors = inputs["CO2 emission intensities"]
    demand_orig_scale = inputs["Demand"]
    # Variability (1s for dispatchable and fraction for VRE)
    availability = inputs["Generation availability"]

    # Storage specific
    if storage_flag
        N = inputs["Storage power to energy ratio"]
        # Efficiency
        F_dch = inputs["Storage discharge loss"]
        F_ch = inputs["Storage charge loss"]
    end

    # Stochasticities
    P = inputs["Demand scenario probabilities"]
    P_f = inputs["Fuel price scenario probabilities"]

    # Representative periods
    n_periods = inputs["Number of representative periods"]
    t_weights = inputs["Period weights"]
    first_periods = inputs["Set of first periods"]
    last_periods = inputs["Set of last periods"]

    # ~~~
    # Build model
    # ~~~

    # Sets
    T = inputs["Number of periods"]
    S = size(P)[1] # number of demand scenarios
    F = size(P_f)[1] # number of fuel cost scenarios
    G = inputs["Number of generation resources"]
    O =  inputs["Number of storage resources"]
    storage_flag = O > 0
    R = G + O


    # ~~~
    # Parameters
    # ~~~
    # Policies
    
    if co2_tax_flag 
        co2_tax = settings["Carbon tax"]
    else
        co2_tax = 0
    end
    if itc_flag
        itc = settings["Tax credits"]
    else
        itc = zeros(R)
    end
    if ptc_flag
        ptc = settings["Production credits"]
    else
        ptc = zeros(R)
    end
    if standard_flag
        standard_perc = settings["Standard mandate"]
    else
        standard_perc = 0
    end

   
    # CVaR parameters
    Ψ = settings["Value-at-Risk percent"] 
    Ω = risk_aversion_weight


    # Scale parameters
    co2_tax = co2_tax/scaling_factor_cost
    demand = demand_orig_scale./scaling_factor_demand
    cost_inv = cost_inv./scaling_factor_cost
    cost_var = cost_var./scaling_factor_cost
    g_max = 50e3/scaling_factor_demand
    # Demand flexibility
    price_nse = [2e3,1800,1100,400]./scaling_factor_cost # originally $/MWh    
    max_nse = [1,0.04,0.024,0.003] # percentage of demand that can be curtailed in each segments
    NSE_SEGS = 4

    # ~~~
    # Model formulation
    # ~~~

    # Generation
    @variable(gep, 0 <= g[r in 1:G, t in 1:T, s in 1:S, f in 1:F] <= g_max) # MWh
    # Capacity
    @variable(gep, 0 <= x[r in 1:R]) # Capacity, MW
    # Non-served energy
    @variable(gep, y[seg in 1:NSE_SEGS, t in 1:T, s in 1:S, f in 1:F] >= 0) # $/MWh
    @expression(gep, eNSE_Cost[t in 1:T, s in 1:S, f in 1:F], sum(price_nse[seg]*y[seg,t,s,f] for seg in 1:NSE_SEGS))

    #if objective_type == "CVaR-Weighted"  ### Always want this so that we can record the CVaR loss and VaR even when not in the objective
        # Auxiliary varliables for CVaR
        @variable(gep, u[s in 1:S, f in 1:F] >= 0) # Gain/loss relative to VaR, $
        @variable(gep, ζ) # VaR variable, $
        @constraint(gep, cvar_tail[s in 1:S, f in 1:F], u[s,f] >= sum(t_weights[t]*g[r,t,s,f]*cost_var[r,f] for r in 1:G, t in 1:T) + sum(t_weights[t]*g[r,t,s,f]*co2_tax*co2_factors[r] for r in 1:G, t in 1:T) + sum(t_weights[t]*eNSE_Cost[t,s,f] for t in 1:T)  - ζ)
    #end 

    if storage_flag
        @variable(gep, z_dch[t in 1:T, s in 1:S, f in 1:F] >= 0)
        @variable(gep, z_ch[t in 1:T, s in 1:S, f in 1:F] >= 0)
        # State of charge
        @variable(gep, e[t in 1:T, s in 1:S, f in 1:F] >= 0) 
        # Primal feasibility
        
        for i in 1:n_periods
            # Wrap first and last periods - start
            @constraints(gep, 
            begin 
                [s in 1:S, f in 1:F], e[first_periods[i],s,f] == e[last_periods[i],s,f] - (1/F_dch)*z_dch[first_periods[i],s,f] + F_ch*z_ch[first_periods[i],s,f]
            end)
            @constraints(gep, 
            begin [s in 1:S, f in 1:F], z_dch[first_periods[i],s,f] <= e[last_periods[i],s,f]
            end)
            # Wrap first and last periods - end
            @constraints(gep, 
            begin [t in first_periods[i]+1:last_periods[i], s in 1:S, f in 1:F], e[t,s,f] == e[t-1,s,f] - (1/F_dch)*z_dch[t,s,f] + F_ch*z_ch[t,s,f] end)
            @constraints(gep, begin [t in first_periods[i]+1:last_periods[i], s in 1:S, f in 1:F], z_dch[t,s,f] <= e[t-1,s,f] end)
        end

        
        # Energy balance for the remaining periods
        @constraint(gep, energy_limit[t in 1:T, s in 1:S, f in 1:F], e[t,s,f] <= (1/N)*x[G+1])
        @constraint(gep, charge_limit_rate[t in 1:T, s in 1:S, f in 1:F], z_ch[t,s,f] <= x[G+1])
        
        @constraint(gep, discharge_limit_rate[t in 1:T, s in 1:S, f in 1:F], z_dch[t,s,f] <= x[G+1])
        @constraint(gep, charge_discharge_balance[t in 1:T, s in 1:S, f in 1:F], z_dch[t,s,f] + z_ch[t,s,f] <= x[G+1])
    end

    # Load shedding
    # Maximum non served energy (primal feasibility)
    @constraint(gep, max_load_shedding[seg in 1:NSE_SEGS,t in 1:T, s in 1:S, f in 1:F], max_nse[seg]*demand[t,s] - y[seg,t,s,f] >= 0)
    @expression(gep,eNSE[t in 1:T,s in 1:S, f in 1:F],sum(y[seg,t,s,f] for seg in 1:NSE_SEGS))
    # Capacity limit on generation
    @expression(gep, capacity_limit_expression[r in 1:G, t in 1:T, s in 1:S, f in 1:F], x[r]*availability[t,r] - g[r,t,s,f])
    @constraint(gep, capacity_limit[r in 1:G, t in 1:T, s in 1:S, f in 1:F], capacity_limit_expression[r,t,s,f] >= 0)

    # Power balance constraint
    if storage_flag
        @constraint(gep, power_balance[t in 1:T, s in 1:S, f in 1:F], sum(g[r,t,s,f] for r in 1:G) + eNSE[t,s,f] + z_dch[t,s,f] - z_ch[t,s,f] - demand[t,s] == 0)
    else
        @constraint(gep, power_balance[t in 1:T, s in 1:S, f in 1:F], sum(g[r,t,s,f] for r in 1:G) + eNSE[t,s,f] - demand[t,s] == 0)
    end

    if standard_flag
        @constraint(gep, standard_constraint[s in 1:S, f in 1:F], sum(t_weights[t]*g[r,t,s,f] for r in 2:3, t in 1:T) >= standard_perc*sum(t_weights[t]*g[r,t,s,f] for r in 1:G, t in 1:T))
    end

    # Define system cost
    @expression(gep, op_cost[s in 1:S, f in 1:F], sum(t_weights[t]*g[r,t,s,f]*cost_var[r,f] for r in 1:G, t in 1:T) + sum(t_weights[t]*eNSE_Cost[t,s,f] for t in 1:T))
    @expression(gep, avg_op_cost[s in 1:S, f in 1:F], op_cost[s,f]/sum(t_weights[t]*demand[t,s] for t in 1:T))
    @expression(gep, inv_cost, sum(x[r]*cost_inv[r] for r in 1:R))
    @expression(gep, avg_inv_cost[s in 1:S], inv_cost/sum(t_weights[t]*demand[t,s] for t in 1:T))
    
    @expression(gep, system_cost, sum(x[r]*cost_inv[r] for r in 1:R) + sum(P[s]*P_f[f]*t_weights[t]*g[r,t,s,f]*cost_var[r,f] for r in 1:G, t in 1:T, s in 1:S, f in 1:F) + sum(P[s]*P_f[f]*t_weights[t]*eNSE_Cost[t,s,f] for t in 1:T, s in 1:S, f in 1:F))

    # Define emissions
    @expression(gep, emissions[s in 1:S, f in 1:F], sum(t_weights[t]*g[r,t,s,f]*co2_factors[r] for r in 1:G, t in 1:T))
    @expression(gep, emissions_exp, sum(P[s]*P_f[f]*emissions[s,f] for s in 1:S, f in 1:F))

    # Generation
    @expression(gep, eANN_GEN[r in 1:G, s in 1:S, f in 1:F], sum(g[r,t,s,f] for t in 1:T))

    
    # Cap constraint
    if cap_flag
        co2_cap = settings["Carbon cap"]
        @constraint(gep, emissions_cap[s in 1:S, f in 1:F], co2_cap >= emissions[s,f])
    else
        co2_cap = 0
    end

    # ~~~
    # Objective function
    # ~~~ 

    # Expectation of system cost
    @expression(gep, sys_cost, sum(x[r]*cost_inv[r]*(1-itc[r]) for r in 1:R) + sum(P[s]*P_f[f]*t_weights[t]*g[r,t,s,f]*(cost_var[r,f]+co2_tax*co2_factors[r]) for r in 1:G, t in 1:T, s in 1:S, f in 1:F) + sum(P[s]*P_f[f]*t_weights[t]*eNSE_Cost[t,s,f] for t in 1:T, s in 1:S, f in 1:F))

    # Weighted CVaR of system cost
    @expression(gep, risk_adjusted_sys_cost, sum(x[r]*cost_inv[r]*(1-itc[r]) for r in 1:R) + Ω*(sum(P[s]*P_f[f]*t_weights[t]*g[r,t,s,f]*(cost_var[r,f]+co2_tax*co2_factors[r]) for r in 1:G, t in 1:T, s in 1:S, f in 1:F) + sum(P[s]*P_f[f]*t_weights[t]*eNSE_Cost[t,s,f] for t in 1:T, s in 1:S, f in 1:F)) + (1-Ω)*(ζ + 1/Ψ*sum(P[s]*P_f[f]*u[s,f] for s in 1:S, f in 1:F)))

    # Sub expression - cvar
    @expression(gep, cvar, (ζ + 1/Ψ*sum(P[s]*P_f[f]*u[s,f] for s in 1:S, f in 1:F)))

    # Sub expression - expected cost
    @expression(gep, expected_cost, sum(P[s]*P_f[f]*t_weights[t]*g[r,t,s,f]*(cost_var[r,f]+co2_tax*co2_factors[r]) for r in 1:G, t in 1:T, s in 1:S, f in 1:F) + sum(P[s]*P_f[f]*t_weights[t]*eNSE_Cost[t,s,f] for t in 1:T, s in 1:S, f in 1:F))

    # Sub expresssion - investment cost
    @expression(gep, investment_cost, sum(x[r]*cost_inv[r]*(1-itc[r]) for r in 1:R))


    if objective_type == "CVaR-Weighted"
        @objective(gep, Min, risk_adjusted_sys_cost)
    else
        # Risk neutral
        @objective(gep, Min, sys_cost)
        
    end
    set_silent(gep)
    gep.ext[:risk_aversion] = obj_weight
    return gep

end

# Run model
function run_optimization_model(gep, inputs, settings)

    optimize!(gep)

    if termination_status(gep) != MOI.OPTIMAL
        @info("Model did not solve to optimality. Status: ", termination_status(gep))
        compute_conflict!(gep)
        list_of_conflicting_constraints = ConstraintRef[];
        for (F, S) in list_of_constraint_types(gep)
            for con in all_constraints(gep, F, S)
                if get_attribute(con, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
                    push!(list_of_conflicting_constraints, con)
                end
            end
        end
        display(list_of_conflicting_constraints)
        error("Model infeasible. See conflicting constraints above.")
    end

     # Numerical settings
    scaling_factor_demand = 1 # to GWs
    scaling_factor_obj = 1 
    scaling_factor_cost = 1 # to $1,000,000s

    # ~~~~
    # Load inputs
    # ~~~

    risk_aversion_flag = settings["Risk aversion flag"]
    risk_aversion_weight = settings["Risk aversion weight"]

    # Settings
    storage_flag = inputs["Storage flag"]
    itc_flag = settings["Investment tax credits flag"]
    ptc_flag = settings["Production tax credits flag"]
    co2_tax_flag = settings["Carbon tax flag"]
    cap_flag = settings["Carbon cap flag"]
    standard_flag = settings["Standard flag"]

    # General
    cost_inv = inputs["Investment costs"]
    cost_var = inputs["Variable costs"]
    co2_factors = inputs["CO2 emission intensities"]
    demand_orig_scale = inputs["Demand"]
    # Variability (1s for dispatchable and fraction for VRE)
    availability = inputs["Generation availability"]

    # Storage specific
    if storage_flag
        N = inputs["Storage power to energy ratio"]
        # Efficiency
        F_dch = inputs["Storage discharge loss"]
        F_ch = inputs["Storage charge loss"]
    end

    # Stochasticities
    P = inputs["Demand scenario probabilities"]
    P_f = inputs["Fuel price scenario probabilities"]

    # Representative periods
    n_periods = inputs["Number of representative periods"]
    t_weights = inputs["Period weights"]
    first_periods = inputs["Set of first periods"]
    last_periods = inputs["Set of last periods"]


    # Sets
    T = inputs["Number of periods"]
    S = size(P)[1] # number of demand scenarios
    F = size(P_f)[1] # number of fuel cost scenarios
    G = inputs["Number of generation resources"]
    O =  inputs["Number of storage resources"]
    storage_flag = O > 0
    R = G + O


    # ~~~
    # Parameters
    # ~~~
    # Policies
    
    if co2_tax_flag 
        co2_tax = settings["Carbon tax"]
    else
        co2_tax = 0
    end
    if itc_flag
        itc = settings["Tax credits"]
    else
        itc = zeros(R)
    end
    if ptc_flag
        ptc = settings["Production credits"]
    else
        ptc = zeros(R)
    end
    if standard_flag
        standard_perc = settings["Standard mandate"]
    else
        standard_perc = 0
    end

    
        # CVaR parameters
        Ψ = settings["Value-at-Risk percent"] 
        Ω = risk_aversion_weight
    

    # Scale parameters
    co2_tax = co2_tax/scaling_factor_cost
    


      # Write outputs
    output = Dict{String, Any}()

        
    price = dual.(gep[:power_balance])
    price_shadow = shadow_price.(gep[:power_balance])  
    mu = dual.(gep[:capacity_limit])
    output["Capacity"] = value.(gep[:x])
    y = value.(gep[:eNSE])
    output["Load shedding"] = value.(gep[:eNSE])
    output["Generation"] = value.(gep[:g])
    output["Capacity dual"] = dual.(gep[:capacity_limit])
    output["Risk adjusted probabilities"] = zeros(S,F)
    output["Consumer surplus"] = zeros(S,F)
    
    theta = dual.(gep[:cvar_tail])/(1-Ω)
    for s in 1:S
        for f in 1:F
            output["Risk adjusted probabilities"][s,f] = theta[s,f]
        end
    end
    


    if storage_flag
        output["Storage charging"] = value.(gep[:z_ch])
        output["Storage discharging"] = value.(gep[:z_dch])
    end
    output["Power balance dual"] = price
    output["Unweighted dual"] = scaling_factor_cost.*price./t_weights
    
    output["Objective function value"] = objective_value(gep)
    output["CVAR duals"] = dual.(gep[:cvar_tail])
    output["System cost"] = value.(gep[:system_cost])
    output["Operating cost"] = value.(gep[:op_cost])
    output["Operating cost average"] = value.(gep[:avg_op_cost])
    output["Investment cost"] = value.(gep[:inv_cost])
    output["Investment cost average"] = value.(gep[:avg_inv_cost])
    output["Scaling factor for costs"] = scaling_factor_cost
    output["Scaling factor for demand"] = scaling_factor_demand
    output["Emissions"] = value.(gep[:emissions])
    output["Emissions expected"] = value.(gep[:emissions_exp])
    # Record results from the linearization or not
    output["Piecewise linearized"] = false
    output["Risk sharing flag"] = true
    output["Optimization model flag"] = true
    
    output["CVaR Loss"] = value.(gep[:u])
    output["VaR"] = value.(gep[:ζ])
    output["CVaR Value"] = value.(gep[:cvar])
    output["Expected cost"] = value.(gep[:expected_cost])
    output["Investment cost"] = value.(gep[:investment_cost])
    output["Risk adjusted system cost"] = value.(gep[:risk_adjusted_sys_cost])
    output["Risk weight"] = gep.ext[:risk_aversion]
    

    if standard_flag
        output["Credit price"] = dual.(gep[:standard_constraint])*1e6
    else
        output["Credit price"] = 0
    end
    if cap_flag
        output["CO2 price"] = dual.(gep[:emissions_cap])*1e6
    else
        output["CO2 price"] = 0
    end

    return output, gep
end

# Set objective

function set_objective!(model, objective_type::String; obj_weight::Float64 = -1.0, set_coeffs::Vector = [])
    if objective_type == "Expectation"
        @objective(model, Min, model[:investment_cost] + model[:expected_cost])
        @info("Objective set to minimize expected system cost: ", objective_function(model))
    elseif objective_type == "Weighted CVaR"
        if obj_weight < -0.0001 || obj_weight > 1.001
            error("For Weighted CVaR objective, please provide a valid weight between 0 and 1.")
        end
        @objective(model, Min, model[:investment_cost] + obj_weight*model[:cvar] + (1-obj_weight)*model[:expected_cost])
        @info("Objective set to minimize weighted CVaR. Weight is: ", obj_weight, " Objective function: ", objective_function(model))
    elseif objective_type == "Capacity"
        if set_coeffs == []
            error("For Capacity objective, please provide a vector of coefficients for the capacity expression.")
        end
        println(set_coeffs)
        @objective(model, Min, set_coeffs'*model[:x])
        @info("Objective set to minimize capacity expression: ", objective_function(model))
    else
        error("Unsupported objective type. Please choose from 'Expectation', 'Weighted CVaR', or 'Capacity'.")
    end
end

function generate_weights(iterations::Int64, num_resources::Int, option::String)
    weights = Vector{Vector{Float64}}(undef, 0)
    vector = zeros(num_resources)
    for i in 1:ceil(Int, iterations/2)
        if option == "random"
            vector = randn(num_resources)*1000
        elseif option == "min/max"
            vector = rand(-1:1, num_resources)*1000
        end
        push!(weights, vector)
        push!(weights, -1*vector) 
    end

    weights = sort_weights(weights)

    return weights
end

function sort_weights(weights::Vector{Vector{Float64}})
    norms = [norm(weights[i]) for i in 1:length(weights)]
    dot_product = collect(dot(weights[i], weights[1]) for i in 1:length(weights))./norms
    dot_product = clamp.(dot_product, -1.0, 1.0)
    angles = [acos(dot_product[i]) for i in 1:length(weights)]
    sorted_indices = sortperm(angles)
    weights = weights[sorted_indices]

    return weights
end

# Add budget constraint

function add_budget_constraint!(model::Model, budget::Float64, budget_type::String; risk_aversion::Float64 = -1.0)
    budget_added = Symbol[]
    if budget_type == "Investment"
        @constraint(model, c_budget_inv, model[:investment_cost] <= budget)
        push!(budget_added, :c_budget_inv)
    elseif budget_type == "Expected"
        @constraint(model, c_budget_exp, model[:expected_cost] <= budget)
        push!(budget_added, :c_budget_exp)
    elseif budget_type == "CVaR"
        @constraint(model, c_budget_cvar, model[:cvar] <= budget)
        push!(budget_added, :c_budget_cvar)
    elseif budget_type == "System_Expected"
        @constraint(model, c_budget_sys_exp, model[:investment_cost] + model[:expected_cost] <= budget)
        push!(budget_added, :c_budget_sys_exp)
    elseif budget_type == "System_CVaR"
        if risk_aversion <= 0 || risk_aversion >= 1
            error("For System_CVaR budget constraint, please provide a valid risk aversion weight between 0 and 1.")
        end
        @constraint(model, c_budget_sys_cvar, model[:investment_cost] + (1-risk_aversion)*model[:cvar] + (risk_aversion)*model[:expected_cost] <= budget)
        push!(budget_added, :c_budget_sys_cvar)
    else
        error("Unsupported budget type. Please choose from 'Investment', 'Expected', 'CVaR', 'System_Expected', or 'System_CVaR'.")
    end

    @info("Added budget constraint: ", model[budget_added[1]])
    return model
end

function add_budget_constraint!(model::Model, budget::Vector, budget_type::Vector; risk_aversion::Float64 = -1.0)
    budget_added = Symbol[]
    for i in 1:length(budget)
        if budget_type[i] == "Investment"
            @constraint(model, c_budget_inv, model[:investment_cost] <= budget[i])
            push!(budget_added, :c_budget_inv)
        elseif budget_type[i] == "Expected"
            @constraint(model, c_budget_exp, model[:expected_cost] <= budget[i])
            push!(budget_added, :c_budget_exp)
        elseif budget_type[i] == "CVaR"
            @constraint(model, c_budget_cvar, model[:cvar] <= budget[i])
            push!(budget_added, :c_budget_cvar)
        elseif budget_type[i] == "System_Expected"
            @constraint(model, c_budget_sys_exp, model[:investment_cost] + model[:expected_cost] <= budget[i])
            push!(budget_added, :c_budget_sys_exp)
        elseif budget_type[i] == "System_CVaR"
            if risk_aversion <= 0 || risk_aversion >= 1
                error("For System_CVaR budget constraint, please provide a valid risk aversion weight between 0 and 1.")
            end
            @constraint(model, c_budget_sys_cvar, model[:investment_cost] + (1-risk_aversion)*model[:cvar] + (risk_aversion)*model[:expected_cost] <= budget[i])
            push!(budget_added, :c_budget_sys_cvar)
        else
            error("Unsupported budget type. Please choose from 'Investment', 'Expected', 'CVaR', 'System_Expected', or 'System_CVaR'.")
        end

        @info("Added budget constraint: ", model[budget_added[i]])
    end
    return model
end

function add_nse_cap!(model::Model, cap::Float64)
    @constraint(model, cNSE_Cap, sum(model[:eNSE][t,s,f] for t in 1:size(model[:eNSE],1), s in 1:size(model[:eNSE],2), f in 1:size(model[:eNSE],3)) <= cap)
end