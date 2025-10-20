# Write results from equilibrium risk-averse GEP model
function write_results(model_output, inputs, settings, results_folder)

    # Parameters
    P = inputs["Demand scenario probabilities"]
    P_f = inputs["Fuel price scenario probabilities"]
    time_index = inputs["Time index"]
    t_weights = inputs["Period weights"]
    co2_factors = inputs["CO2 emission intensities"]
    scaling_factor_cost= model_output["Scaling factor for costs"]
    scaling_factor_demand = model_output["Scaling factor for demand"]

    # Variables
    cap = model_output["Capacity"] 
    shed = model_output["Load shedding"]
    gen = model_output["Generation"]
    m = model_output["Capacity dual"]
    theta = model_output["Risk adjusted probabilities"]
    Ω = settings["Risk aversion weight"]
    price = model_output["Power balance dual"]
    
    # Settings
    resources = inputs["Generation resources"]
    all_resources = inputs["Resources"]
    storage_flag = inputs["Storage flag"]

    # Sets
    T = inputs["Number of periods"]
    S = size(P)[1] # number of demand scenarios
    F = size(P_f)[1] # number of fuel cost scenarios
    G = inputs["Number of generation resources"]
    O =  inputs["Number of storage resources"]
    R = G + O 


    if storage_flag
        charge = model_output["Storage charging"]
        discharge = model_output["Storage discharging"]
        net_battery = charge - discharge
    end

    # Collect results
    df_cap = DataFrame(Resource = all_resources, Capacity = cap)
    # VaR
    if !model_output["Risk sharing flag"]
        value_at_risk = model_output["VaR"] # VaR is a vector indexed by technology
    else
        value_at_risk = ones().*model_output["VaR"] # VaR is a single value here
    end
    df_var = DataFrame(Resource = all_resources, VaR = value_at_risk)
    # Weights for periods
    df_weights = DataFrame()
    insertcols!(df_weights, 1, :Time => time_index)
    insertcols!(df_weights, :Weights => t_weights)
    # Collect time series results per scenario
    df_bat = DataFrame()
    insertcols!(df_bat, 1, :Time => time_index)
    df_gen = DataFrame()
    insertcols!(df_gen, 1, :Time => time_index)
    df_price = DataFrame()
    insertcols!(df_price, 1, :Time => time_index)
    
    df_co2_all = DataFrame()
    df_nse = DataFrame()
    insertcols!(df_nse, 1, :Time => time_index)
    df_u = DataFrame(Resource = all_resources)
    df_theta = DataFrame()
    insertcols!(df_theta, 1, :Resource => all_resources)
    df_oper_all = DataFrame()

    # Total system cost in x$
    system_cost = model_output["System cost"]
    # Operating cost expected in x$
    exp_op_cost = sum(P[s]*P_f[f]*model_output["Operating cost"][s,f] for s in 1:S, f in 1:F)
    # Average operating cost expected in x$/MWh
    exp_avg_op_cost = sum(P[s]*P_f[f]*model_output["Operating cost average"][s,f] for s in 1:S, f in 1:F)
    # Average invesment cost expected in x$/mWh
    exp_avg_inv_cost = sum(P[s]*P_f[f]*model_output["Investment cost average"][s] for s in 1:S, f in 1:F)
    # Average system cost in x$/MWh
    avg_syscost = exp_avg_op_cost + exp_avg_inv_cost
    # Save to dataframe
    df_syscost = DataFrame(SystemCost = system_cost, AverageCost=avg_syscost)
    # df_avg_syscost = DataFrame(SystemCost = avg_syscost)
    
    # Risk-adjusted system cost
    # For now only works when a single scenario is in the CVaR tail
    cvar = maximum(model_output["Operating cost"])
    avg_cvar = maximum(model_output["Operating cost average"])
    sys_cost_risk = model_output["Investment cost"] + Ω*exp_op_cost + (1-Ω)*cvar
    avg_sys_cost_risk = model_output["Investment cost"] + Ω*exp_avg_op_cost + (1-Ω)*avg_cvar
    # Save to dataframe
    df_syscost_risk = DataFrame(SystemCostRisk = sys_cost_risk, AverageCostRisk = avg_sys_cost_risk)

    # Objective 
    objective_val = model_output["Objective function value"]
    df_obj = DataFrame(Objective = objective_val)

    # Emissions
    co2 = model_output["Emissions"]
    co2_exp = model_output["Emissions expected"]
    df_emissions = DataFrame(Expected_emissions = co2_exp)

    # Collect scenario results
    for s in 1:S
        for f in 1:F
            # Generation
            col_names = [string(i,"_Demand-",string(s),"_FuelPrice-",string(f)) for i in resources]
            df_gen = hcat(df_gen, DataFrame(transpose(gen[:,:,s,f]), col_names))
            # Power price
            col_name = string("Demand-",string(s),"_FuelPrice-",string(f))
            insertcols!(df_price, col_name => price[:,s,f])

            # Collect CO2 emissions
            insertcols!(df_co2_all, col_name => model_output["Emissions"][s,f])
            # Collect operating cost
            insertcols!(df_oper_all, col_name => model_output["Operating cost"][s,f])
            
            if !model_output["Risk sharing flag"]
                insertcols!(df_theta, col_name => theta[:,s,f])
            else
                insertcols!(df_theta, col_name => theta[s,f])
            end
            
            # Battery operation
            if storage_flag
                insertcols!(df_bat, col_name => net_battery[:,s,f])
            end
            # Load shedding
            insertcols!(df_nse, col_name => shed[:,s,f])
            
            if !model_output["Risk sharing flag"]
                insertcols!(df_u, col_name => model_output["CVaR Loss"][:,s,f])
            else
                insertcols!(df_u, col_name => model_output["CVaR Loss"][s,f])
            end
        end
    end

    if Sys.isunix()
        sep = "/"
    elseif Sys.iswindows()
        sep = "\U005c"
    end

    # Write output files
    CSV.write(string(results_folder,sep,"capacity.csv"), df_cap)
    CSV.write(string(results_folder,sep,"generation.csv"), df_gen)
    if storage_flag
        CSV.write(string(results_folder,sep,"battery.csv"), df_bat)
    end

    CSV.write(string(results_folder,sep,"price.csv"), df_price)
    CSV.write(string(results_folder,sep,"weights.csv"), df_weights)
    CSV.write(string(results_folder,sep,"nse.csv"), df_nse)
    CSV.write(string(results_folder,sep,"objective.csv"), df_obj)
    CSV.write(string(results_folder,sep,"emissions.csv"), df_emissions)
    CSV.write(string(results_folder,sep,"emissions_all.csv"), df_co2_all)
    CSV.write(string(results_folder,sep,"oper_cost_all.csv"), df_oper_all)
    if settings["Risk aversion flag"] == true
        CSV.write(string(results_folder,sep,"VaR.csv"), df_var)
        CSV.write(string(results_folder,sep,"loss_cvar.csv"), df_u)
        CSV.write(string(results_folder,sep,"risk-adj_probs.csv"), df_theta)
    end
    CSV.write(string(results_folder,sep,"system_cost.csv"), df_syscost)
    CSV.write(string(results_folder,sep,"risk_system_cost.csv"), df_syscost_risk)

end


