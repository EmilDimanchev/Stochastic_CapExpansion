# Write results from equilibrium risk-averse GEP model
function write_results(model_output, inputs, settings, results_folder)

    if !isdir(results_folder)
        mkpath(results_folder)
    end

    # Parameters
    P = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
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
    # Average system cost by scenario
    avg_sys_cost_scenarios = [model_output["Operating cost average"][s,f] + model_output["Investment cost average"][s] for s in 1:S, f in 1:F]
    # Save to dataframe
    df_syscost = DataFrame(SystemCost = system_cost, AverageCost=avg_syscost)
    # df_avg_syscost = DataFrame(SystemCost = avg_syscost)
    
    # Risk-adjusted system cost
    # For now only works when a single scenario is in the CVaR tail
    #cvar = maximum(model_output["Operating cost"])
    avg_cvar = maximum(model_output["Operating cost average"])
    sys_cost_risk = model_output["Risk adjusted system cost"] #model_output["Investment cost"] + Ω*exp_op_cost + (1-Ω)*cvar
    cvar_op = maximum(model_output["Operating cost"])
    sys_cost_risk_expost = model_output["Investment cost"] + Ω*exp_op_cost + (1-Ω)*cvar_op
    avg_sys_cost_risk = model_output["Investment cost"] + Ω*exp_avg_op_cost + (1-Ω)*avg_cvar
    println("Investment cost: ", model_output["Investment cost"])
    println("Expected operating cost: ", exp_avg_op_cost)
    println("Average CVaR operating cost: ", avg_cvar)
    println("Model CVaR: ", model_output["CVaR Value"])
    println("Worst-case operating cost: ", cvar_op)
    # Save to dataframe
    df_syscost_risk = DataFrame(SystemCostRisk = sys_cost_risk, SystemCostRisk_ExPost = sys_cost_risk_expost, AverageCostRisk = avg_sys_cost_risk)

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

            # Average system cost by scenario
            col_name_syscost = string("Average_System_Cost_", "Demand-",string(s),"_FuelPrice-",string(f))
            insertcols!(df_syscost, col_name_syscost => avg_sys_cost_scenarios[s,f])

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

    return df_cap, df_syscost, df_syscost_risk, df_emissions, df_oper_all

end

function write_results_benders(results::Dict, inputs::Dict, settings::Dict, results_folder::String; budgets = nothing)

    if !isdir(results_folder)
    mkpath(results_folder)
    end
    # This function is for writing results from the Benders algorithm, which has a different structure than the monolithic model output
    MP_output = results["MP"]
    SPs_output = results["SPs"]

    # Should only have all uncertainties set to false if SPs eval are provided
    eval_SPs = []
    try 
        eval_SPs = (haskey(results, "SPs_eval")) ? results["SPs_eval"] : []
    catch
        error("SPs_eval not found in results. Make sure to include SPs_eval in the output of the Benders algorithm when running with all uncertainties set to false.")
    end 

    if !isdir(results_folder)
        mkpath(results_folder)
    end

    # Parameters
    P_s = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]
    P_s_eval = [0.0]
    P_f_eval = [0.0]
    P_k_eval = [0.0]
    

    if haskey(results, "SPs_eval")
        P_s_eval = inputs["Full Demand scenario probabilities"]
        P_f_eval = inputs["Full Gas price scenario probabilities"]
        P_k_eval = inputs["Full Weather scenario probabilities"]
    end
    S_eval = size(P_s_eval)[1] # number of demand scenarios
    F_eval = size(P_f_eval)[1] # number of fuel cost scenarios
    K_eval = size(P_k_eval)[1] # number of weather scenarios

    time_index = inputs["Time index"]
    t_weights = inputs["Period weights"]
    co2_factors = inputs["CO2 emission intensities"]
    scaling_factor_cost= settings["Scaling factor cost"]
    scaling_factor_demand = settings["Scaling factor demand"]

    # Variables
    cap = MP_output["Capacity"]
 
    shed = nothing
    gen = nothing
    
    #m = model_output["Capacity dual"]
    #theta = model_output["Risk adjusted probabilities"]
    Ω = settings["Risk aversion weight"]
    price = nothing
    # Settings
    resources = inputs["Generation resources"]
    all_resources = inputs["Resources"]
    storage_flag = inputs["Storage flag"]

    # Sets
    T = inputs["Number of periods"]
    S = size(P_s)[1] # number of demand scenarios
    F = size(P_f)[1] # number of fuel cost scenarios
    K = size(P_k)[1] # number of weather scenarios
    G = inputs["Number of generation resources"]
    O =  inputs["Number of storage resources"]
    R = G + O 
    L = inputs["Number of lines"]
    Z = inputs["Number of zones"]


    #if storage_flag
     #   charge = collect(SPs_output[s,f,k]["Storage charging"] for s in 1:size(P_s)[1], f in 1:size(P_f)[1], k in 1:size(P_k)[1])
     #   discharge = collect(SPs_output[s,f,k]["Storage discharging"] for s in 1:size(P_s)[1], f in 1:size(P_f)[1], k in 1:size(P_k)[1])
     #   net_battery = charge - discharge
    #end

    # Collect results
    df_cap = DataFrame(Resource = all_resources, Capacity = cap)
    # VaR
    if settings["Risk sharing flag"]
        value_at_risk = MP_output["VaR"] # VaR is a vector indexed by technology
    else
        value_at_risk = ones().*MP_output["VaR"] # VaR is a single value here
    end
    df_var = DataFrame(Resource = all_resources, VaR = value_at_risk)
    # Weights for periods
    df_weights = DataFrame()
    insertcols!(df_weights, 1, :Time => time_index)
    insertcols!(df_weights, :Weights => reshape(t_weights, :))
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

    df_syscost = DataFrame()
    df_emissions = DataFrame()
    inv_cost_by_zone = MP_output["Inv cost by zone"]

    if eval_SPs != [] #### Condition means all uncertainties off and we have Eval SPs
        # Base expressions
        system_cost_ev_base = MP_output["Inv_cost"] + results["Expected Value"] # model_output["System cost"] Base system cost for single scenario
        system_cost_cvar_base = MP_output["Inv_cost"] + Ω*results["Expected Value"] + (1-Ω)*results["CVaR"] # model_output["Risk adjusted system cost"] Base risk-adjusted system cost for single scenario
        cvar_op_base = results["CVaR"] # model_output["CVaR Value"] Base CVaR for single scenario
        ev_op_base = results["Expected Value"] # model_output["Expected operating cost"] Base expected operating cost for single scenario
        co2 = collect(SPs_output[s,f,k]["Emissions"] for s in 1:size(P_s)[1], f in 1:size(P_f)[1], k in 1:size(P_k)[1]) #model_output["Emissions"]
        co2_exp = sum(P_s[s]*P_f[f]*P_k[k]*co2[s,f,k] for s in 1:S, f in 1:F, k in 1:K)
        prices = collect(SPs_output[s,f,k]["Power price"] for s in 1:size(P_s)[1], f in 1:size(P_f)[1], k in 1:size(P_k)[1])
        prices_exp = sum(P_s[s]*P_f[f]*P_k[k]*sum(prices[s,f,k]) for s in 1:S, f in 1:F, k in 1:K)

        #by zone
        cost_by_zone_base = collect(SPs_output[s,f,k]["Cost by zone"] for s in 1:size(P_s)[1], f in 1:size(P_f)[1], k in 1:size(P_k)[1]) 
        CO2_by_zone_base = collect(SPs_output[s,f,k]["Emissions by zone"] for s in 1:size(P_s)[1], f in 1:size(P_f)[1], k in 1:size(P_k)[1])
        cost_by_zone_exp = sum(P_s[s]*P_f[f]*P_k[k]*cost_by_zone_base[s,f,k] for s in 1:S, f in 1:F, k in 1:K) .+ inv_cost_by_zone
        co2_by_zone_exp = sum(P_s[s]*P_f[f]*P_k[k]*CO2_by_zone_base[s,f,k] for s in 1:S, f in 1:F, k in 1:K)
        prices_by_zone_exp = sum(P_s[s]*P_f[f]*P_k[k]*sum(prices[s,f,K][t,:] for t in 1:T) for s in 1:size(P_s)[1], f in 1:size(P_f)[1], k in 1:size(P_k)[1])/length(prices[1,1,1][:,1])

        # Eval SP expressions
         # Total system cost in x$
        system_cost_eval = results["Expected Value Eval"] + MP_output["Inv_cost"]
        # Operating cost expected in x$
        exp_op_cost_eval = results["Expected Value Eval"]
        sys_cost_risk_eval = MP_output["Inv_cost"] + Ω*results["Expected Value Eval"] + (1-Ω)*results["CVaR_Eval"] #model_output["Investment cost"] + Ω*exp_op_cost + (1-Ω)*cvar
        cvar_op_eval = results["CVaR_Eval"]
        co2_eval = collect(eval_SPs[s,f,k]["Emissions"] for s in 1:size(P_s_eval)[1], f in 1:size(P_f_eval)[1], k in 1:size(P_k_eval)[1]) #model_output["Emissions"]
        co2_exp_eval = sum(P_s_eval[s]*P_f_eval[f]*P_k_eval[k]*co2_eval[s,f,k] for s in 1:S, f in 1:F, k in 1:K) #model_output["Emissions expected"]
        prices_eval = collect(eval_SPs[s,f,k]["Power price"] for s in 1:size(P_s_eval)[1], f in 1:size(P_f_eval)[1], k in 1:size(P_k_eval)[1])
        prices_exp_eval = sum(P_s_eval[s]*P_f_eval[f]*P_k_eval[k]*sum(prices_eval[s,f,k]) for s in 1:S, f in 1:F, k in 1:K)/length(prices_eval[1,1,1]) #model_output["Expected Power price"]

        #by zone
        cost_by_zone_eval = collect(eval_SPs[s,f,k]["Cost by zone"] for s in 1:size(P_s_eval)[1], f in 1:size(P_f_eval)[1], k in 1:size(P_k_eval)[1])
        CO2_by_zone_eval = collect(eval_SPs[s,f,k]["Emissions by zone"] for s in 1:size(P_s_eval)[1], f in 1:size(P_f_eval)[1], k in 1:size(P_k_eval)[1])
        cost_by_zone_eval_exp = sum(P_s_eval[s]*P_f_eval[f]*P_k_eval[k]*cost_by_zone_eval[s,f,k] for s in 1:S, f in 1:F, k in 1:K) .+ inv_cost_by_zone
        co2_by_zone_eval_exp = sum(P_s_eval[s]*P_f_eval[f]*P_k_eval[k]*CO2_by_zone_eval[s,f,k] for s in 1:S, f in 1:F, k in 1:K)
        prices_by_zone_eval_exp = sum(P_s_eval[s]*P_f_eval[f]*P_k_eval[k]*sum(prices_eval[s,f,k][t,:] for t in 1:T) for s in 1:S, f in 1:F, k in 1:K)/length(prices_eval[1,1,1][:,1])
        
        # Make dfs
        df_syscost = DataFrame(Investment_Cost = MP_output["Inv_cost"], System_Cost_EV_Base = system_cost_ev_base, System_Cost_CVaR_Base = system_cost_cvar_base, CVaR_OpCost_Base = cvar_op_base, EV_OpCost_Base = ev_op_base, System_Cost_EV_Eval = system_cost_eval, System_Cost_CVaR_Eval = sys_cost_risk_eval, CVaR_OpCost_Eval = cvar_op_eval, EV_OpCost_Eval = exp_op_cost_eval, Expected_Power_Price_Eval = prices_exp_eval, Expected_Power_Price_Base = prices_exp)        
        df_emissions = DataFrame(Expected_CO2_Base = co2_exp, Expected_CO2_Eval = co2_exp_eval)
        for z in 1:Z
            col_name_cost_zone = string("Cost_by_zone_base", string("Zone"), string(z))
            insertcols!(df_syscost, col_name_cost_zone => cost_by_zone_exp[z])
            col_name_co2_zone = string("CO2_by_zone_base", string("Zone"), string(z))
            insertcols!(df_emissions, col_name_co2_zone => co2_by_zone_exp[z])
            col_name_price_zone = string("Price_by_zone_base", string("Zone"), string(z))
            insertcols!(df_syscost, col_name_price_zone => prices_by_zone_exp[z])
            
            col_name_cost_zone = string("Cost_by_zone_eval", string("Zone"), string(z))
            insertcols!(df_syscost, col_name_cost_zone => cost_by_zone_eval_exp[z])
            col_name_co2_zone = string("CO2_by_zone_eval", string("Zone"), string(z))
            insertcols!(df_emissions, col_name_co2_zone => co2_by_zone_eval_exp[z])
            col_name_price_zone = string("Price_by_zone_eval", string("Zone"), string(z))
            insertcols!(df_syscost, col_name_price_zone => prices_by_zone_eval_exp[z])
        end
        if settings["Write all scenarios flag"]
            for s in 1:size(P_s_eval)[1]
                for f in 1:size(P_f_eval)[1]
                    for k in 1:size(P_k_eval)[1]
                        # Collect CO2 emissions
                        col_name_Emissions = string("EmissionsD",string(s),"F",string(f),"W",string(k))
                        insertcols!(df_co2_all, col_name_Emissions => eval_SPs[s,f,k]["Emissions"])
                        # Collect operating cost
                        col_name_OpCost = string("OperatingCostD",string(s),"F",string(f),"W",string(k))
                        insertcols!(df_syscost, col_name_OpCost => eval_SPs[s,f,k]["SP objective"])

                        # Collect power prices
                        col_name_PowerPrice = string("avgPowerPriceD",string(s),"F",string(f),"W",string(k))
                        insertcols!(df_syscost, col_name_PowerPrice => sum(eval_SPs[s,f,k]["Power price"])/length(eval_SPs[1,1,1]["Power price"]))

                    end
                end
            end
        end

    else #### Proceed as normal
         # Total system cost in x$
        system_cost = results["Expected Value"] + MP_output["Inv_cost"]
        # Operating cost expected in x$
        exp_op_cost = results["Expected Value"]
        sys_cost_risk = MP_output["Inv_cost"] + Ω*exp_op_cost + (1-Ω)*results["CVaR"] #model_output["Investment cost"] + Ω*exp_op_cost + (1-Ω)*cvar
        cvar_op = results["CVaR"]
        prices = collect(SPs_output[s,f,k]["Power price"] for s in 1:size(P_s)[1], f in 1:size(P_f)[1], k in 1:size(P_k)[1])
        prices_exp = sum(P_s[s]*P_f[f]*P_k[k]*sum(prices[s,f,k]) for s in 1:S, f in 1:F, k in 1:K)/length(prices[1,1,1]) #model_output["Expected Power price"]
        

        cost_by_zone_base = collect(SPs_output[s,f,k]["Cost by zone"] for s in 1:size(P_s)[1], f in 1:size(P_f)[1], k in 1:size(P_k)[1])
        CO2_by_zone_base = collect(SPs_output[s,f,k]["Emissions by zone"] for s in 1:size(P_s)[1], f in 1:size(P_f)[1], k in 1:size(P_k)[1])
        cost_by_zone_exp = sum(P_s[s]*P_f[f]*P_k[k]*cost_by_zone_base[s,f,k] for s in 1:S, f in 1:F, k in 1:K)
        co2_by_zone_exp = sum(P_s[s]*P_f[f]*P_k[k]*CO2_by_zone_base[s,f,k] for s in 1:S, f in 1:F, k in 1:K)
        prices_by_zone_exp = sum(P_s[s]*P_f[f]*P_k[k]*sum(prices[s,f,K][t,:] for t in 1:T) for s in 1:size(P_s)[1], f in 1:size(P_f)[1], k in 1:size(P_k)[1])/length(prices[1,1,1][:,1])

        
       
        co2 = collect(SPs_output[s,f,k]["Emissions"] for s in 1:size(P_s)[1], f in 1:size(P_f)[1], k in 1:size(P_k)[1]) #model_output["Emissions"]
        co2_exp = sum(P_s[s]*P_f[f]*P_k[k]*co2[s,f,k] for s in 1:S, f in 1:F, k in 1:K) #model_output["Emissions expected"]
        df_syscost = DataFrame(Investment_Cost = MP_output["Inv_cost"], EV_SystemCost = system_cost, CVaR_SystemCost= sys_cost_risk, CVaR_OpCost = cvar_op, EV_OpCost = exp_op_cost, Expected_Power_Price = prices_exp)
        df_emissions = DataFrame(Expected_emissions = co2_exp)
        for z in 1:Z
            col_name_cost_zone = string("Cost_by_zone_base", string("Zone"), string(z))
            insertcols!(df_syscost, col_name_cost_zone => cost_by_zone_exp[z])
            col_name_co2_zone = string("CO2_by_zone_base", string("Zone"), string(z))
            insertcols!(df_emissions, col_name_co2_zone => co2_by_zone_exp[z])
            col_name_price_zone = string("Price_by_zone_base", string("Zone"), string(z))
            insertcols!(df_syscost, col_name_price_zone => prices_by_zone_exp[z])
        end
        if settings["Write all scenarios flag"]
            for s in 1:size(P_s)[1]
                for f in 1:size(P_f)[1]
                    for k in 1:size(P_k)[1]
                        # Collect CO2 emissions
                        col_name_Emissions = string("EmissionsD",string(s),"F",string(f),"W",string(k))
                        insertcols!(df_co2_all, col_name_Emissions => SPs_output[s,f,k]["Emissions"])
                        # Collect operating cost
                        col_name_OpCost = string("OperatingCostD",string(s),"F",string(f),"W",string(k))
                        insertcols!(df_syscost, col_name_OpCost => SPs_output[s,f,k]["SP objective"])
                        # Collect power prices
                        col_name_PowerPrice = string("avgPowerPriceD",string(s),"F",string(f),"W",string(k))
                        insertcols!(df_syscost, col_name_PowerPrice => sum(SPs_output[s,f,k]["Power price"])/length(SPs_output[1,1,1]["Power price"]))
                    end
                end
            end
        end
    end
   

    if settings["Write all scenarios flag"]
        shed = eval_SPs == [] ? collect(SPs_output[s,f,k]["Load shedding"] for s in 1:size(P_s)[1], f in 1:size(P_f)[1], k in 1:size(P_k)[1]) : collect(eval_SPs[s,f,k]["Load shedding"] for s in 1:size(P_s_eval)[1], f in 1:size(P_f_eval)[1], k in 1:size(P_k_eval)[1])
        gen = eval_SPs == [] ? collect(SPs_output[s,f,k]["Generation"] for s in 1:size(P_s)[1], f in 1:size(P_f)[1], k in 1:size(P_k)[1]) : collect(eval_SPs[s,f,k]["Generation"] for s in 1:size(P_s_eval)[1], f in 1:size(P_f_eval)[1], k in 1:size(P_k_eval)[1])
        price = eval_SPs == [] ? collect(SPs_output[s,f,k]["Power price"] for s in 1:size(P_s)[1], f in 1:size(P_f)[1], k in 1:size(P_k)[1]) : collect(eval_SPs[s,f,k]["Power price"] for s in 1:size(P_s_eval)[1], f in 1:size(P_f_eval)[1], k in 1:size(P_k_eval)[1])
        # Collect scenario results
        for s in 1:S
            for f in 1:F
                for k in 1:K
                    # Generation
                    col_names = [string(i,"_D",string(s),"F",string(f),"W",string(k)) for i in resources]
                    df_gen = hcat(df_gen, DataFrame(transpose(gen[s,f,k]), col_names))
                    for z in 1:Z
                        
                        # Power price
                        col_name = string("D",string(s),"F",string(f),"W",string(k), "Z", string(z))
                        insertcols!(df_price, col_name => price[s,f,k][:,z])

                        # Average system cost by scenario
                        #col_name_syscost = string("Average_System_Cost_", "Demand-",string(s),"_FuelPrice-",string(f),"_Weather-",string(k))
                        #insertcols!(df_syscost, col_name_syscost => avg_sys_cost_scenarios[s,f,k])

                        # Collect CO2 emissions
                        #insertcols!(df_co2_all, col_name => model_output["Emissions"][s,f,k])
                        # Collect operating cost
                        #insertcols!(df_oper_all, col_name => model_output["Operating cost"][s,f,k])
                        
                        #if !model_output["Risk sharing flag"]
                        #   insertcols!(df_theta, col_name => theta[:,s,f,k])
                        #else
                        #    insertcols!(df_theta, col_name => theta[s,f,k])
                        #end
                        
                        # Battery operation
                        #if storage_flag
                        #    insertcols!(df_bat, col_name => net_battery[:,s,f,k])
                        #end
                        # Load shedding
                        insertcols!(df_nse, col_name => shed[s,f,k][:,z])
                        
                        #if !model_output["Risk sharing flag"]
                        #    insertcols!(df_u, col_name => model_output["CVaR Loss"][:,s,f,k])
                        #else
                        #    insertcols!(df_u, col_name => model_output["CVaR Loss"][s,f,k])
                        #end
                    end
                end
            end
        end
    end

    if Sys.isunix()
        sep = "/"
    elseif Sys.iswindows()
        sep = "\U005c"
    end

    names = inputs["Resources"]
    gaps = results["Gaps"]
    if (typeof(gaps[1]) == Vector{Float64}) && !isnothing(budgets)
        k = collect(key for key in keys(budgets))
        names = vcat(k, names)
    else
        names = ["Gap"; names]
    end
    CSV.write(joinpath(results_folder, "Convergence.csv"), DataFrame(hcat(stack(gaps, dims=1), stack(results["Capacity per iteration"], dims=1)), names))

    if settings["Write all scenarios flag"]
        # Write output files
        #CSV.write(string(results_folder,sep,"capacity.csv"), df_cap)
        CSV.write(string(results_folder,sep,"generation.csv"), df_gen)
        #if storage_flag
        #    CSV.write(string(results_folder,sep,"battery.csv"), df_bat)
        #end

        CSV.write(string(results_folder,sep,"price.csv"), df_price)
       # CSV.write(string(results_folder,sep,"weights.csv"), df_weights)
        CSV.write(string(results_folder,sep,"nse.csv"), df_nse)
        #CSV.write(string(results_folder,sep,"objective.csv"), df_obj)
        #CSV.write(string(results_folder,sep,"emissions.csv"), df_emissions)
        #CSV.write(string(results_folder,sep,"emissions_all.csv"), df_co2_all)
        #CSV.write(string(results_folder,sep,"oper_cost_all.csv"), df_oper_all)
        #if settings["Risk aversion flag"] == true
          #  CSV.write(string(results_folder,sep,"VaR.csv"), df_var)
          #  CSV.write(string(results_folder,sep,"loss_cvar.csv"), df_u)
          #  CSV.write(string(results_folder,sep,"risk-adj_probs.csv"), df_theta)
        #end
        #CSV.write(string(results_folder,sep,"system_cost.csv"), df_syscost)
        #CSV.write(string(results_folder,sep,"risk_system_cost.csv"), df_syscost_risk)

    end

    return df_cap, df_syscost, df_emissions   

end

function write_gaps!(gaps::AbstractVector, labels::AbstractVector,results_folder::String)
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    df_gaps = DataFrame(gaps, labels)
    CSV.write(joinpath(results_folder, "Gaps.csv"), df_gaps)
end


function write_exploration_results!(dfs_cap::AbstractVector, dfs_syscost::AbstractVector, dfs_emissions::AbstractVector, results_folder::String, labels::AbstractVector = [], scenario::String = "")

    if !isdir(results_folder)
        mkpath(results_folder)
    end

    names = gather_names(dfs_cap[1], dfs_syscost[1], dfs_emissions[1], length(labels) > 0)
    data = Array{Any,2}(undef, length(dfs_cap), length(names))
    data = gather_data(data, dfs_cap, dfs_syscost, dfs_emissions, labels)
    df_exploration = DataFrame(data, names)

    CSV.write(joinpath(results_folder,"Summary_"*scenario*".csv"), df_exploration)

end

function gather_names(cap::DataFrame, sys_cost::DataFrame, emissions::DataFrame, label_bool::Bool = false)
    name_vec = []
    for resource in cap.Resource
        push!(name_vec, resource)
    end
    for col in names(sys_cost)
        push!(name_vec, col)
    end
    for col in names(emissions)
        push!(name_vec, col)
    end
    if label_bool
        name_vec = vcat("Iteration Name", name_vec)
    end
    return name_vec
end

function gather_data(data::AbstractArray,dfs_cap::AbstractVector, dfs_syscost::AbstractVector, dfs_emissions::AbstractVector, labels::AbstractVector)
    
    for i in 1:length(dfs_cap)
        row = make_row(dfs_cap[i], dfs_syscost[i], dfs_emissions[i], length(labels) > 0 ? labels[i] : "")
        data[i,:] = row
    end

    return data
end

function make_row(cap::DataFrame, sys_cost::DataFrame, emissions::DataFrame, label::String)
    row = Any[]
    if label != ""
        push!(row, label)
    end
    for resource in cap.Resource
        push!(row, cap[cap.Resource .== resource, :Capacity][1])
    end
    for col in names(sys_cost)
        push!(row, sys_cost[1,col])
    end
    for col in names(emissions)
        push!(row, emissions[1,col])
    end

    return row
end