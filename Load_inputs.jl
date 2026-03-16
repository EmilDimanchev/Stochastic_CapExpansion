function load_input_data(inputs_path, settings)

    # Collect data here
    inputs = Dict{String, Any}()

    if Sys.isunix()
        sep = "/"
    elseif Sys.iswindows()
        sep = "\U005c"
    end

    # Settings
    demand_risk = settings["Demand risk flag"]
    fuel_risk = settings["Fuel risk flag"]

    # Input parameters
    period_weights_input = CSV.read(string(inputs_path,sep,"Representative_period_weights.csv"), DataFrame, header=true)
    # Technologies
    resources_input = CSV.read(string(inputs_path,sep,"Resources.csv"), DataFrame, header=true)
    
    inputs["Set of CES resources"] = findall(x->x==1, resources_input[!, "CES_eligible"])
    inputs["Set of non-CES resources"] = findall(x->x==0, resources_input[!, "CES_eligible"])

    gen_resources = resources_input[1:end,1]
    inputs["Generation resources"] = gen_resources
    x_max = resources_input[1:end,"Capacity_UB"]
    resource_names = copy(gen_resources)
    
    storage_input = CSV.read(string(inputs_path,sep,"Resources_storage.csv"), DataFrame, header=true)
    storage_flag = size(storage_input)[1] > 0
    var_cost_input = CSV.read(string(inputs_path,sep,"Variable_cost.csv"), DataFrame, header=true)
    resource_avail_input = CSV.read(string(inputs_path,sep,"Resources_availability.csv"), DataFrame, header=true)
    # # Piecewise linearization inputs
    # linearization_input = CSV.read(string(inputs_path,sep,"Linearization_segments.csv"), DataFrame, header=true)
    # Stochasticities
    probabilities_input_demand = CSV.read(string(inputs_path,sep,"Scenario_probabilities_demand.csv"), DataFrame, header=true)
    probabilities_input_fuel = CSV.read(string(inputs_path,sep,"Scenario_probabilities_fuel.csv"), DataFrame, header=true)


    # Demand
    if demand_risk
        demand_input = CSV.read(string(inputs_path,sep,"Demand.csv"), DataFrame, header=true)
        inputs["Demand"] = round.(Array(demand_input[:,3:end]),digits=1)
        inputs["Demand scenario probabilities"] = Array(probabilities_input_demand[probabilities_input_demand[!,:Uncertainty] .== "Demand",2:end][1,:])
    else
        demand_input = CSV.read(string(inputs_path,sep,"Demand.csv"), DataFrame, header=true, select=["Time_index", "Demand-base"])
        inputs["Demand"] = round.(Array(demand_input[:,2]),digits=1)
        inputs["Demand scenario probabilities"] = [1]
    end
    time_index = demand_input[:,1]
    inputs["Time index"] = time_index

    # Fuel price
    if fuel_risk
        # Variable cost in $/MWh
        inputs["Variable costs"] = Array(var_cost_input[1:end, 3:end])
        inputs["Fuel price scenario probabilities"] = Array(probabilities_input_fuel[probabilities_input_fuel[!,:Uncertainty] .== "Fuel",2:end][1,:])
    else
        inputs["Variable costs"] = Array(var_cost_input[1:end, 2]) # $/MWh
        inputs["Fuel price scenario probabilities"] = [1]
    end

    inputs["Number of generation resources"] = size(resources_input)[1]
    inputs["Number of storage resources"] = size(storage_input)[1]
    
    # Storage
    if storage_flag
        # inputs["Storage power to energy ratio"] = storage_input[1:end,"Power_to_energy_ratio"][1]
        # inputs["Storage discharge loss"] = storage_input[1:end,"Loss_discharge"][1]
        # inputs["Storage charge loss"] = storage_input[1:end,"Loss_charge"][1]
        # max_stor_cap = storage_input[1:end,"Capacity_UB"][1]
        # push!(x_max, max_stor_cap) 
        # push!(resource_names, "Battery") 
        inputs["Storage power to energy ratio"] = storage_input[1:end,"Power_to_energy_ratio"][1]
        inputs["Storage discharge loss"] = storage_input[1:end,"Loss_discharge"][1]
        inputs["Storage charge loss"] = storage_input[1:end,"Loss_charge"][1]
        stor_resources = storage_input[1:end,1]
        max_stor = storage_input[1:end,"Capacity_UB"]
        x_max_all = vcat(x_max, max_stor)
        resources_all = vcat(resource_names, stor_resources)
        inputs["Resources"] = resources_all
    else
        inputs["Resources"] = resource_names
        x_max_all = x_max
    end   

    availability = round.(Matrix(resource_avail_input[:, 2:end]),digits=2)
    inputs["Generation availability"] = availability
    inputs["Capacity upper bounds"] = x_max_all
    # inputs["Piecewise linearization inputs"] = linearization_input

    # Investment cost
    cost_inv = round.(resources_input[1:end,"Investment cost"],digits=2) # $/MW-year
    # Add storage cost
    if storage_flag
        cost_inv = vcat(cost_inv, round.(storage_input[1:end,"Investment cost"],digits=2)) 
    end
    co2_factors = round.(resources_input[:,"Emissions_ton_per_MWh"],digits=2) # ton per MWh

    # Temporal domain
    n_periods = size(period_weights_input)[1]
    inputs["Number of representative periods"] = n_periods
    T_length = size(demand_input)[1]
    inputs["Number of periods"] = T_length

    # Representative period weights
    t_weights = zeros(T_length)
    period_length = Integer(T_length/n_periods)
    first_periods = collect(1:period_length:T_length) 
    last_periods = collect(period_length:period_length:T_length)
    inputs["Set of first periods"] = first_periods
    inputs["Set of last periods"] = last_periods

    for i in 1:n_periods
        t_weights[first_periods[i]:last_periods[i]] .= period_weights_input[!,"Weight"][i]/period_length
    end

    inputs["Period weights"] = t_weights
    inputs["Investment costs"] = cost_inv
    inputs["CO2 emission intensities"] = co2_factors
    inputs["Storage flag"] = storage_flag
    inputs["Storage"] = storage_input

    return inputs
end

