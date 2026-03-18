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
    price_weighted = model_output["Power balance dual"]
    price = scaling_factor_cost.*price_weighted./t_weights
    
    # Settings
    resources = inputs["Generation resources"]
    all_resources = inputs["Resources"]
    linearized = model_output["Piecewise linearized"]
    optimization_model_flag = model_output["Optimization model flag"]
    storage_flag = inputs["Storage flag"]

    # Sets
    T = inputs["Number of periods"]
    S = size(P)[1] # number of demand scenarios
    F = size(P_f)[1] # number of fuel cost scenarios
    G = inputs["Number of generation resources"]
    O =  inputs["Number of storage resources"]
    R = G + O 

    # Revenues
    marginal_revenues_gen = model_output["Generation marginal revenues"]
    if storage_flag
        marginal_revenues_stor = model_output["Storage marginal revenues"]
    end
    marginal_revenues = zeros(R,S,F)
    total_revenues = zeros(R,S,F)
    for r in 1:R
        for s in 1:S
            for f in 1:F
                if r <= G
                    marginal_revenues[r,s,f] = marginal_revenues_gen[r,s,f]
                    total_revenues[r,s,f] = marginal_revenues_gen[r,s,f]*cap[r]
                else
                    total_revenues[r,s,f] = marginal_revenues_stor[s,f]*cap[r]
                    marginal_revenues[r,s,f] = marginal_revenues_stor[s,f]
                end
            end
        end
    end

    if linearized
        linearized_rev = model_output["Linearized total revenues"]
        slack = model_output["Linearization slack"]
        revenue_error_abs = linearized_rev.-total_revenues
        revenue_error_prc = replace!(replace!(replace!(abs.(revenue_error_abs)./total_revenues, NaN=>0)), Inf=>0)
        revenue_error_max = maximum(revenue_error_prc)

        revenue_error = replace(linearized_rev./total_revenues .- 1, -Inf=>0)

        inv_error = (Ω*sum(P[s]*P_f[f]*slack[:,s,f] for s in 1:S, f in 1:F) + (1 - Ω).*sum(theta[:,s,f].*slack[:,s,f] for s in 1:S, f in 1:F))./(Ω*sum(P[s]*P_f[f]*total_revenues[:,s,f] for s in 1:S, f in 1:F) + (1 - Ω).*sum(theta[:,s,f,].*total_revenues[:,s,f] for s in 1:S, f in 1:F))
        iso_error = sum(slack[r,:,:] for r in 1:R)./sum(total_revenues[r,:,:] for r in 1:R)

        df_slacks = DataFrame()
        insertcols!(df_slacks, 1, :Resource => all_resources)
        df_slacks_perc = DataFrame()
        insertcols!(df_slacks_perc, 1, :Resource => all_resources)
        df_inv_error = DataFrame(Resource = all_resources, Error = inv_error)
    end

    if storage_flag
        charge = model_output["Storage charging"]
        discharge = model_output["Storage discharging"]
        net_battery = charge - discharge

        if !optimization_model_flag
            soc_dual_weighted = model_output["SOC dual"]
            soc_dual = scaling_factor_cost.*soc_dual_weighted./t_weights
            ch_dual = scaling_factor_cost.*model_output["Charging dual"]./t_weights
            dch_dual = scaling_factor_cost.*model_output["Discharging dual"]./t_weights
            dch_e_dual = scaling_factor_cost.*model_output["Discharging energy dual"]./t_weights
            bal_dual = scaling_factor_cost.*model_output["Balance dual"]./t_weights
            energy_dual = scaling_factor_cost.*model_output["Energy capacity dual"]./t_weights
            
            df_soc_dual = DataFrame()
            insertcols!(df_soc_dual, 1, :Time => time_index)
            df_ch_dual = DataFrame()
            insertcols!(df_ch_dual, 1, :Time => time_index)
            df_dch_dual = DataFrame()
            insertcols!(df_dch_dual, 1, :Time => time_index)
            df_dch_e_dual = DataFrame()
            insertcols!(df_dch_e_dual, 1, :Time => time_index)
            df_bal_dual = DataFrame()
            insertcols!(df_bal_dual, 1, :Time => time_index)
            df_energy_dual = DataFrame()
            insertcols!(df_energy_dual, 1, :Time => time_index)
        end

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
    df_revenues = DataFrame(Resource = all_resources)
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
    # co2_exp = sum(P[s]*P_f[f]*co2 for s in 1:S, f in 1:F) 
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
            if !optimization_model_flag
                if storage_flag
                    insertcols!(df_soc_dual, col_name => soc_dual[:,s,f])
                    insertcols!(df_ch_dual, col_name => ch_dual[:,s,f])
                    insertcols!(df_dch_dual, col_name => dch_dual[:,s,f])
                    insertcols!(df_dch_e_dual, col_name => dch_e_dual[:,s,f])
                    insertcols!(df_bal_dual, col_name => bal_dual[:,s,f])
                    insertcols!(df_energy_dual, col_name => energy_dual[:,s,f])
                end
            end

            # Collect CO2 emissions
            insertcols!(df_co2_all, col_name => model_output["Emissions"][s,f])
            # Collect operating cost
            insertcols!(df_oper_all, col_name => model_output["Operating cost"][s,f])
            
            if !model_output["Risk sharing flag"]
                insertcols!(df_theta, col_name => theta[:,s,f])
            else
                insertcols!(df_theta, col_name => theta[s,f])
            end

            if linearized
                insertcols!(df_slacks, col_name => slack[:,s,f])
                insertcols!(df_slacks_perc, col_name => revenue_error_prc[:,s,f])
            end
            # Emissions
            
            # Battery operation
            if storage_flag
                insertcols!(df_bat, col_name => net_battery[:,s,f])
            end
            # Load shedding
            insertcols!(df_nse, col_name => shed[:,s,f])
            insertcols!(df_revenues, col_name => marginal_revenues[:,s,f])
            if !model_output["Risk sharing flag"]
                insertcols!(df_u, col_name => model_output["CVaR Loss"][:,s,f])
            else
                insertcols!(df_u, col_name => model_output["CVaR Loss"][s,f])
            end
        end
    end


    if !(isdir(string(results_folder)))
        mkdir(string(results_folder))
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
    if linearized
        CSV.write(string(results_folder,sep,"investment_slack.csv"), df_inv_error)
        CSV.write(string(results_folder,sep,"slack_absolute.csv"), df_slacks)
        CSV.write(string(results_folder,sep,"slack_perc.csv"), df_slacks_perc)
    end
    CSV.write(string(results_folder,sep,"price.csv"), df_price)
    if storage_flag
        if !optimization_model_flag
            CSV.write(string(results_folder,sep,"phi_soc.csv"), df_soc_dual)
            CSV.write(string(results_folder,sep,"phi_ch.csv"), df_ch_dual)
            CSV.write(string(results_folder,sep,"phi_dch.csv"), df_dch_dual)
            CSV.write(string(results_folder,sep,"xi_d.csv"), df_dch_e_dual)
            CSV.write(string(results_folder,sep,"phi_bal.csv"), df_bal_dual)
            CSV.write(string(results_folder,sep,"phi_cap.csv"), df_energy_dual)
        end
    end
    CSV.write(string(results_folder,sep,"weights.csv"), df_weights)
    CSV.write(string(results_folder,sep,"nse.csv"), df_nse)
    CSV.write(string(results_folder,sep,"objective.csv"), df_obj)
    CSV.write(string(results_folder,sep,"emissions.csv"), df_emissions)
    CSV.write(string(results_folder,sep,"emissions_all.csv"), df_co2_all)
    CSV.write(string(results_folder,sep,"oper_cost_all.csv"), df_oper_all)
    CSV.write(string(results_folder,sep,"revenues.csv"), df_revenues)
    if settings["Risk aversion flag"] == true
        CSV.write(string(results_folder,sep,"VaR.csv"), df_var)
        CSV.write(string(results_folder,sep,"loss_cvar.csv"), df_u)
        CSV.write(string(results_folder,sep,"risk-adj_probs.csv"), df_theta)
    end
    CSV.write(string(results_folder,sep,"system_cost.csv"), df_syscost)
    CSV.write(string(results_folder,sep,"risk_system_cost.csv"), df_syscost_risk)

end


