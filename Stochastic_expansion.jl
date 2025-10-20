using JuMP
using CSV
using DataFrames
using Gurobi

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
    # price_cap = settings["Price cap"]/scaling_factor_cost # originally $/MWh before scaling
    price_nse = settings["Price cap"]./scaling_factor_cost # originally $/MWh    
    max_nse = [1] # percentage of demand that can be curtailed in each segments
    NSE_SEGS = 1
    # price_nse = [2e3,1800,1100,400]./scaling_factor_cost # originally $/MWh    
    # max_nse = [1,0.04,0.024,0.003] # percentage of demand that can be curtailed in each segments
    # NSE_SEGS = 4


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
        @variable(gep, u[s in 1:S, f in 1:F] >= 0) # loss relative to VaR, $/MW
        @variable(gep, ζ) # VaR variable, $/MW
        # @constraint(gep, cvar_tail[s in 1:S, f in 1:F], u[s,f] >= sum(t_weights[t]*g[r,t,s,f]*cost_var[r,f] for r in 1:G, t in 1:T) + sum(t_weights[t]*eNSE_Cost[t,s,f] for t in 1:T) + sum(cost_inv[r]*x[r] for r in 1:R) - ζ)
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

    # if settings["Emission test"]
    #     @constraint(gep, emish_test, sys_cost <= 2.819466116)
    #     # @constraint(gep, emish_test, 2.819466116 >= sys_cost)
    # end

    if risk_aversion_flag
        @objective(gep, Min, sum(x[r]*cost_inv[r]*(1-itc[r]) for r in 1:R) + Ω*(sum(P[s]*P_f[f]*t_weights[t]*g[r,t,s,f]*(cost_var[r,f]+co2_tax*co2_factors[r]) for r in 1:G, t in 1:T, s in 1:S, f in 1:F) + sum(P[s]*P_f[f]*t_weights[t]*eNSE_Cost[t,s,f] for t in 1:T, s in 1:S, f in 1:F)) + (1-Ω)*(ζ + 1/Ψ*sum(P[s]*P_f[f]*u[s,f] for s in 1:S, f in 1:F)))

        # eps = 0.000001
        # @constraint(gep, mga, sum(x[r]*cost_inv[r]*(1-itc[r]) for r in 1:R) + Ω*(sum(P[s]*P_f[f]*t_weights[t]*g[r,t,s,f]*(cost_var[r,f]+co2_tax*co2_factors[r]) for r in 1:G, t in 1:T, s in 1:S, f in 1:F) + sum(P[s]*P_f[f]*t_weights[t]*eNSE_Cost[t,s,f] for t in 1:T, s in 1:S, f in 1:F)) + (1-Ω)*(ζ + 1/Ψ*sum(P[s]*P_f[f]*u[s,f] for s in 1:S, f in 1:F)) <= 3.447079778 + eps)

        # @constraint(gep, mga_lb, sum(x[r]*cost_inv[r]*(1-itc[r]) for r in 1:R) + Ω*(sum(P[s]*P_f[f]*t_weights[t]*g[r,t,s,f]*(cost_var[r,f]+co2_tax*co2_factors[r]) for r in 1:G, t in 1:T, s in 1:S, f in 1:F) + sum(P[s]*P_f[f]*t_weights[t]*eNSE_Cost[t,s,f] for t in 1:T, s in 1:S, f in 1:F)) + (1-Ω)*(ζ + 1/Ψ*sum(P[s]*P_f[f]*u[s,f] for s in 1:S, f in 1:F)) >= 3.447079778-eps)

        # @objective(gep, Max, emissions_exp)
    else
        @objective(gep, Min, sys_cost)
        
    end

    optimize!(gep)

    # Write outputs
    output = Dict{String, Any}()

        
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
                output["Consumer surplus"][s,f] = sum(t_weights.*price_nse.*demand[:,s] - price[:,s,f].*demand[:,s] + y[:,s,f].*(price[:,s,f] - t_weights.*price_nse))
            end
        end
    end

    # output["--- risk-adjusted CS"] = output["Consumer surplus"][1,1]*0.75+output["Consumer surplus"][2,1]*0.25

    # Calculate revenues
    revenues = zeros(R,S,F)
    for r in 1:R
        for s in 1:S
            for f in 1:F
                if r <= G
                    # revenues[r,s,f] = sum(mu[r,t,s,f]*availability[t,r] for t in 1:T)/(Ω*P[s]*P_f[f]+(1-Ω)*theta[s,f])
                    revenues[r,s,f] = sum(mu[r,t,s,f]*availability[t,r] for t in 1:T)
                    
                    
                    # revenues[r,s,f] = sum((price[t,s,f]/(Ω*P[s]*P_f[f]+(1-Ω)*theta[s,f])-cost_var[r,f])*value.(g)[r,t,s,f]/value.(x)[r] for t in 1:T)
                    
                else
                    # Ignore storage revenues for now
                end
            end
        end
    end

    
    output["Generation marginal revenues"] = revenues
    output["Storage marginal revenues"] = zeros(S,F)
    if storage_flag
        output["Storage charging"] = value.(z_ch)
        output["Storage discharging"] = value.(z_dch)
    end
    output["Power balance dual"] = price
    output["Unweighted dual"] = scaling_factor_cost.*price./t_weights
    
    # if settings["Emission test"]
    #     power_price = settings["power price 1"]
    #     cost_dual = dual.(emish_test)
    #     cost_shadow_dual = shadow_price.(emish_test)
    #     output["new_dual"] = cost_dual
    #     output["new_shadow"] = cost_shadow_dual
    #     # output["Unweighted dual_2"] = (price .- power_price.*cost_dual)./t_weights
    #     output["Marginal emissions"] = (price .- power_price.*cost_dual)./t_weights
    #     output["Power balance dual consistent"] = (price_shadow .- power_price.*cost_shadow_dual)./t_weights
    # end
    
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