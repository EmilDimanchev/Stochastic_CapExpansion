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

#=====================================================================

Shared results schema for the Benders / MGA / mapping pipeline

`make_results_dfs` is the single source of truth for how a solved
capacity vector (whether it comes from a live MP/SP solve, an
interior-sampled point, or an interpolated point) gets turned into the
three summary dataframes (`df_cap`, `df_syscost`, `df_emissions`) that
feed both the per-run CSVs (`write_results_benders`) and the
cross-run summary (`write_exploration_results!`). Everything downstream
depends on these three dataframes always sharing the same columns
regardless of whether `mapping` is true or false, or whether
out-of-sample evaluation subproblems (`eval_SPs`) are supplied.

======================================================================#

function make_results_dfs(cap::AbstractVector, line_cap::AbstractVector, tot_inv_cost, inv_cost_by_zone::AbstractVector,
                           SPs_output, cvar, ev, inputs::Dict, settings::Dict;
                           eval_SPs = [], cvar_eval = nothing, ev_eval = nothing)

    # Parameters
    P_s = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]
    Ω = settings["Risk aversion weight"]

    # Sets
    T = inputs["Number of periods"]
    S = size(P_s)[1] # number of demand scenarios
    F = size(P_f)[1] # number of fuel cost scenarios
    K = size(P_k)[1] # number of weather scenarios
    L = inputs["Number of lines"]
    Z = inputs["Number of zones"]
    line_names = collect("Line " .* string.(1:L))
    all_names = vcat(inputs["Resources"], line_names)

    # ----- Capacity (wide format: resources + lines as columns) -----
    all_cap = vcat(cap, line_cap)
    df_cap = DataFrame(all_cap', all_names)

    # ----- Base (in-sample) system cost / emissions / prices -----
    system_cost = ev + tot_inv_cost
    sys_cost_risk = tot_inv_cost + Ω*ev + (1-Ω)*cvar

    prices = collect(SPs_output[s,f,k]["Power price"] for s in 1:S, f in 1:F, k in 1:K)
    prices_exp = sum(P_s[s]*P_f[f]*P_k[k]*sum(prices[s,f,k]) for s in 1:S, f in 1:F, k in 1:K)/length(prices[1,1,1])

    cost_by_zone_base = collect(SPs_output[s,f,k]["Cost by zone"] for s in 1:S, f in 1:F, k in 1:K)
    CO2_by_zone_base = collect(SPs_output[s,f,k]["Emissions by zone"] for s in 1:S, f in 1:F, k in 1:K)
    cost_by_zone_exp = sum(P_s[s]*P_f[f]*P_k[k]*cost_by_zone_base[s,f,k] for s in 1:S, f in 1:F, k in 1:K)
    co2_by_zone_exp = sum(P_s[s]*P_f[f]*P_k[k]*CO2_by_zone_base[s,f,k] for s in 1:S, f in 1:F, k in 1:K)
    prices_by_zone_exp = sum(P_s[s]*P_f[f]*P_k[k]*sum(prices[s,f,k][t,:] for t in 1:T) for s in 1:S, f in 1:F, k in 1:K)/length(prices[1,1,1][:,1])

    co2 = collect(SPs_output[s,f,k]["Emissions"] for s in 1:S, f in 1:F, k in 1:K)
    co2_exp = sum(P_s[s]*P_f[f]*P_k[k]*co2[s,f,k] for s in 1:S, f in 1:F, k in 1:K)

    df_syscost = DataFrame(Investment_Cost = tot_inv_cost, EV_SystemCost = system_cost, CVaR_SystemCost = sys_cost_risk,
                            CVaR_OpCost = cvar, EV_OpCost = ev, Expected_Power_Price = prices_exp)
    df_emissions = DataFrame(Expected_emissions = co2_exp)

    for z in 1:Z
        insertcols!(df_syscost, string("Cost_by_zone_base", "Zone", z) => cost_by_zone_exp[z])
        insertcols!(df_emissions, string("CO2_by_zone_base", "Zone", z) => co2_by_zone_exp[z])
        insertcols!(df_syscost, string("Price_by_zone_base", "Zone", z) => prices_by_zone_exp[z])
    end

    if settings["Write all scenarios flag"]
        for s in 1:S, f in 1:F, k in 1:K
            insertcols!(df_emissions, string("Emissions_D_", s, "_F_", f, "_W_", k) => SPs_output[s,f,k]["Emissions"])
            insertcols!(df_syscost, string("OperatingCost_D_", s, "_F_", f, "_W_", k) => SPs_output[s,f,k]["SP objective"])
            for z in 1:Z
                insertcols!(df_syscost, string("avgPowerPrice_by_zone_D_", s, "_F_", f, "_W_", k, "_Zone", z) => sum(SPs_output[s,f,k]["Power price"][:,z])/length(SPs_output[1,1,1]["Power price"][:,z]))
            end
        end
    end

    # ----- Out-of-sample eval columns, only when eval_SPs is supplied -----
    if eval_SPs != []
        cvar_eval === nothing && error("cvar_eval must be provided when eval_SPs is non-empty")
        ev_eval === nothing && error("ev_eval must be provided when eval_SPs is non-empty")

        P_s_eval = inputs["Full Demand scenario probabilities"]
        P_f_eval = inputs["Full Gas price scenario probabilities"]
        P_k_eval = inputs["Full Weather scenario probabilities"]
        S_eval = size(P_s_eval)[1]
        F_eval = size(P_f_eval)[1]
        K_eval = size(P_k_eval)[1]

        system_cost_eval = ev_eval + tot_inv_cost
        sys_cost_risk_eval = tot_inv_cost + Ω*ev_eval + (1-Ω)*cvar_eval

        prices_eval = collect(eval_SPs[s,f,k]["Power price"] for s in 1:S_eval, f in 1:F_eval, k in 1:K_eval)
        prices_exp_eval = sum(P_s_eval[s]*P_f_eval[f]*P_k_eval[k]*sum(prices_eval[s,f,k]) for s in 1:S_eval, f in 1:F_eval, k in 1:K_eval)/length(prices_eval[1,1,1])

        cost_by_zone_eval = collect(eval_SPs[s,f,k]["Cost by zone"] for s in 1:S_eval, f in 1:F_eval, k in 1:K_eval)
        CO2_by_zone_eval = collect(eval_SPs[s,f,k]["Emissions by zone"] for s in 1:S_eval, f in 1:F_eval, k in 1:K_eval)
        cost_by_zone_eval_exp = sum(P_s_eval[s]*P_f_eval[f]*P_k_eval[k]*cost_by_zone_eval[s,f,k] for s in 1:S_eval, f in 1:F_eval, k in 1:K_eval)
        co2_by_zone_eval_exp = sum(P_s_eval[s]*P_f_eval[f]*P_k_eval[k]*CO2_by_zone_eval[s,f,k] for s in 1:S_eval, f in 1:F_eval, k in 1:K_eval)
        prices_by_zone_eval_exp = sum(P_s_eval[s]*P_f_eval[f]*P_k_eval[k]*sum(prices_eval[s,f,k][t,:] for t in 1:T) for s in 1:S_eval, f in 1:F_eval, k in 1:K_eval)/length(prices_eval[1,1,1][:,1])

        co2_eval = collect(eval_SPs[s,f,k]["Emissions"] for s in 1:S_eval, f in 1:F_eval, k in 1:K_eval)
        co2_exp_eval = sum(P_s_eval[s]*P_f_eval[f]*P_k_eval[k]*co2_eval[s,f,k] for s in 1:S_eval, f in 1:F_eval, k in 1:K_eval)

        insertcols!(df_syscost, :EV_SystemCost_Eval => system_cost_eval)
        insertcols!(df_syscost, :CVaR_SystemCost_Eval => sys_cost_risk_eval)
        insertcols!(df_syscost, :CVaR_OpCost_Eval => cvar_eval)
        insertcols!(df_syscost, :EV_OpCost_Eval => ev_eval)
        insertcols!(df_syscost, :Expected_Power_Price_Eval => prices_exp_eval)
        insertcols!(df_emissions, :Expected_emissions_Eval => co2_exp_eval)

        for z in 1:Z
            insertcols!(df_syscost, string("Cost_by_zone_eval", "Zone", z) => cost_by_zone_eval_exp[z])
            insertcols!(df_emissions, string("CO2_by_zone_eval", "Zone", z) => co2_by_zone_eval_exp[z])
            insertcols!(df_syscost, string("Price_by_zone_eval", "Zone", z) => prices_by_zone_eval_exp[z])
        end

        if settings["Write all scenarios flag"]
            for s in 1:S_eval, f in 1:F_eval, k in 1:K_eval
                insertcols!(df_emissions, string("Emissions_Eval_D_", s, "_F_", f, "_W_", k) => eval_SPs[s,f,k]["Emissions"])
                insertcols!(df_syscost, string("OperatingCost_Eval_D_", s, "_F_", f, "_W_", k) => eval_SPs[s,f,k]["SP objective"])
                for z in 1:Z
                    insertcols!(df_syscost, string("avgPowerPrice_by_zone_Eval_D_", s, "_F_", f, "_W_", k, "_Zone", z) => sum(eval_SPs[s,f,k]["Power price"][:,z])/length(eval_SPs[1,1,1]["Power price"][:,z]))
                end
            end
        end
    end

    return df_cap, df_syscost, df_emissions
end

function make_results_mapping_dfs(MP_output::Dict, SPs_output, cvar, ev, inputs::Dict, settings::Dict;
                                   eval_SPs = [], cvar_eval = nothing, ev_eval = nothing)
    return make_results_dfs(MP_output["Capacity"], MP_output["Line expansion"], MP_output["Inv_cost"], MP_output["Inv cost by zone"],
                             SPs_output, cvar, ev, inputs, settings; eval_SPs = eval_SPs, cvar_eval = cvar_eval, ev_eval = ev_eval)
end

function make_results_mapping_dfs(capacities::AbstractVector, line_capacities::AbstractVector, SPs_output, tot_inv_cost, inv_cost_zone, cvar, ev, inputs::Dict, settings::Dict;
                                   eval_SPs = [], cvar_eval = nothing, ev_eval = nothing)
    return make_results_dfs(capacities, line_capacities, tot_inv_cost, inv_cost_zone,
                             SPs_output, cvar, ev, inputs, settings; eval_SPs = eval_SPs, cvar_eval = cvar_eval, ev_eval = ev_eval)
end

function write_results_benders(results::Dict, inputs::Dict, settings::Dict, results_folder::String; budgets = nothing)

    if !isdir(results_folder)
        mkpath(results_folder)
    end

    MP_output = results["MP"]
    SPs_output = results["SPs"]
    eval_SPs = haskey(results, "SPs_eval") ? results["SPs_eval"] : []
    cvar_eval = haskey(results, "CVaR_Eval") ? results["CVaR_Eval"] : nothing
    ev_eval = haskey(results, "Expected Value Eval") ? results["Expected Value Eval"] : nothing

    # Standard summary dataframes, identical in structure to the mapping=true pipeline
    df_cap, df_syscost, df_emissions = make_results_mapping_dfs(MP_output, SPs_output, results["CVaR"], results["Expected Value"], inputs, settings;
                                                                  eval_SPs = eval_SPs, cvar_eval = cvar_eval, ev_eval = ev_eval)

    # Parameters needed for the raw scenario time-series output below
    P_s = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]
    time_index = inputs["Time index"]
    resources = inputs["Generation resources"]

    S = size(P_s)[1] # number of demand scenarios
    F = size(P_f)[1] # number of fuel cost scenarios
    K = size(P_k)[1] # number of weather scenarios
    Z = inputs["Number of zones"]

    if eval_SPs != []
        P_s_eval = inputs["Full Demand scenario probabilities"]
        P_f_eval = inputs["Full Gas price scenario probabilities"]
        P_k_eval = inputs["Full Weather scenario probabilities"]
        S_use, F_use, K_use = size(P_s_eval)[1], size(P_f_eval)[1], size(P_k_eval)[1]
        scenario_source = eval_SPs
    else
        S_use, F_use, K_use = S, F, K
        scenario_source = SPs_output
    end

    if Sys.isunix()
        sep = "/"
    elseif Sys.iswindows()
        sep = "\U005c"
    end

    names = inputs["Resources"]
    gaps = results["Gaps"]
    if (typeof(gaps[1]) == Vector{Float64}) && !isnothing(budgets)
        budget_keys = collect(key for key in keys(budgets))
        names = vcat(budget_keys, names)
    else
        names = ["Gap"; names]
    end
    names = haskey(results, "Time Reg hist") ? [names; "MP Solve Time"; "SP Solve Time"; "Reg Solve Time"] : [names; "MP Solve Time"; "SP Solve Time"]
    data = haskey(results, "Time Reg hist") ? hcat(stack(gaps, dims=1), stack(results["Capacity per iteration"], dims=1), stack(results["Time MP hist"], dims=1), stack(results["Time SP hist"], dims=1), stack(results["Time Reg hist"], dims=1)) : hcat(stack(gaps, dims=1), stack(results["Capacity per iteration"], dims=1), stack(results["Time MP hist"], dims=1), stack(results["Time SP hist"], dims=1))
    CSV.write(joinpath(results_folder, "Convergence.csv"), DataFrame(data, names))

    if settings["Write all scenarios flag"]
        df_gen = DataFrame()
        insertcols!(df_gen, 1, :Time => time_index)
        df_price = DataFrame()
        insertcols!(df_price, 1, :Time => time_index)
        df_nse = DataFrame()
        insertcols!(df_nse, 1, :Time => time_index)

        for s in 1:S_use, f in 1:F_use, k in 1:K_use
            col_names = [string(i,"_D",s,"F",f,"W",k) for i in resources]
            df_gen = hcat(df_gen, DataFrame(transpose(scenario_source[s,f,k]["Generation"]), col_names))
            for z in 1:Z
                col_name = string("D",s,"F",f,"W",k,"Z",z)
                insertcols!(df_price, col_name => scenario_source[s,f,k]["Power price"][:,z])
                insertcols!(df_nse, col_name => scenario_source[s,f,k]["Load shedding"][:,z])
            end
        end

        CSV.write(string(results_folder,sep,"generation.csv"), df_gen)
        CSV.write(string(results_folder,sep,"price.csv"), df_price)
        CSV.write(string(results_folder,sep,"nse.csv"), df_nse)
    end

    return df_cap, df_syscost, df_emissions

end

function remove_undefs(x)
    x_new = []
    for i in eachindex(x)
        if !isassigned(x, i)
            break
        end
        push!(x_new, x[i])
    end
    return x_new
end

function write_mapping_results(results::Dict, inputs, settings)
    mp_by_iteration = remove_undefs(results["All MP outputs per iteration"])
    sp_by_iteration = remove_undefs(results["All SP outputs per iteration"])
    cvars = remove_undefs(results["CVaR_hist"])
    evs = remove_undefs(results["Expected Value hist"])

    n_iterations = length(mp_by_iteration)

    mp_output = mp_by_iteration[1]
    sp_output = sp_by_iteration[1]
    cvar = cvars[1]
    ev = evs[1]

    df_cap, df_syscost, df_emissions = make_results_mapping_dfs(mp_output, sp_output, cvar, ev, inputs, settings)

    for i in 2:n_iterations
        mp_output = mp_by_iteration[i]
        sp_output = sp_by_iteration[i]
        cvar = cvars[i]
        ev = evs[i]
        temp_df_cap, temp_df_syscost, temp_df_emissions = make_results_mapping_dfs(mp_output, sp_output, cvar, ev, inputs, settings)
        append!(df_cap, temp_df_cap)
        append!(df_syscost, temp_df_syscost)
        append!(df_emissions, temp_df_emissions)
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

function write_exploration_results!(results_cap::AbstractVector, results_syscost::AbstractVector, results_emissions::AbstractVector, results_folder::String, run_labels::AbstractVector = [], summary_name::String = ""; mapping = false)

    if !isdir(results_folder)
        mkpath(results_folder)
    end

    # df_cap/df_syscost/df_emissions share the same wide-format schema regardless of
    # whether they came from write_results_benders or make_results_mapping_dfs (mapping=true
    # or mapping=false), so a single vcat + join works for every caller.
    dfs_cap = vcat(results_cap...)
    dfs_syscost = vcat(results_syscost...)
    dfs_emissions = vcat(results_emissions...)

    total_length = nrow(dfs_cap)
    # run_labels only has a 1:1 correspondence with the rows when every pushed result was
    # given a label (e.g. mapping=false, or a fully-labeled interpolate evaluation). When
    # mapping=true adds unlabeled interior samples on top of the labeled runs, fall back to
    # generic sample labels.
    labels = length(run_labels) == total_length ? String.(run_labels) : ["Sample "*string(i) for i in 1:total_length]

    dfs_cap.Index = labels
    dfs_syscost.Index = labels
    dfs_emissions.Index = labels

    df_exploration = innerjoin(dfs_cap, dfs_syscost, on = :Index)
    df_exploration = innerjoin(df_exploration, dfs_emissions, on = :Index)

    CSV.write(joinpath(results_folder,"Summary_"*summary_name*".csv"), df_exploration)

end
