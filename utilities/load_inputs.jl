function load_input_data(inputs_path, settings)

    # Collect data here
    inputs = Dict{String, Any}()


    # Technologies
    resources_input = CSV.read(string(inputs_path,"/","Resources.csv"), DataFrame, header=true)

    gen_resources = resources_input[1:end,1]
    inputs["Generation resources"] = gen_resources
    x_max = resources_input[1:end,"Capacity_UB"]
    resource_names = copy(gen_resources)
    
    storage_input = CSV.read(string(inputs_path,"/","Resources_storage.csv"), DataFrame, header=true)
    storage_flag = size(storage_input)[1] > 0
    var_cost_input = CSV.read(string(inputs_path,"/","Variable_cost.csv"), DataFrame, header=true)
    

    # Keep track of number of resources
    inputs["Number of generation resources"] = size(resources_input)[1]
    inputs["Number of storage resources"] = size(storage_input)[1]

    # Stochasticities
    probabilities_input_demand = CSV.read(string(inputs_path,"/","Scenario_probabilities_demand.csv"), DataFrame, header=true)
    probabilities_input_fuel = CSV.read(string(inputs_path,"/","Scenario_probabilities_fuel.csv"), DataFrame, header=true)
    probabilities_input_weather = CSV.read(string(inputs_path,"/","Scenario_probabilities_weather.csv"), DataFrame, header=true)

    # Demand
    demand_input = CSV.read(string(inputs_path,"/","Demand.csv"), DataFrame, header=true)
    
    time_index = demand_input[:,1]
    inputs["Time index"] = time_index
    demand_add_input = CSV.read(string(inputs_path,"/","Demand_adders.csv"), DataFrame, header=true)

    if settings["Demand uncertainty"]
        inputs["Demand scenario probabilities"] = Array(probabilities_input_demand[probabilities_input_demand[!,:Uncertainty] .== "Demand",2:end][1,:])
        demand_add = Array(demand_add_input[1,2:end])
    else
        inputs["Demand scenario probabilities"] = [1]
        demand_add = [0]
    end
    
    inputs["Demand adders"] = demand_add
    # Central demand scenario
    inputs["Demand"] = round.(Array(demand_input[:,2]),digits=1)
    # Deterministic case uses expected demand
    inputs["Demand shift weather"] = zeros(length(inputs["Demand"]))

    if settings["Weather uncertainty"]
        
        inputs["Demand shift weather"] = round.(Array(demand_input[:,3:end]) .- Array(demand_input[:,2]), digits=1)
        
        inputs["Weather scenario probabilities"] = Array(probabilities_input_weather[probabilities_input_weather[!,:Uncertainty] .== "Weather",2:end][1,:])
        # Availability
        
        availability = zeros(length(inputs["Time index"]), inputs["Number of generation resources"], length(inputs["Weather scenario probabilities"]))
        
        for scen in 1:length(inputs["Weather scenario probabilities"])
            resource_avail_input = CSV.read(string(inputs_path,"/","Resources_availability_",string(scen),".csv"), DataFrame, header=true)
            availability[:,:,scen] = Array(resource_avail_input[:, 2:end])
        end
    else

        inputs["Weather scenario probabilities"] = [1]
        # Availability
        resource_avail_input = CSV.read(string(inputs_path,"/","Resources_availability_1.csv"), DataFrame, header=true)
        availability = round.(Matrix(resource_avail_input[:, 2:end]),digits=2)

        inputs["Weather scenario probabilities"] = [1]
    end
    
    inputs["Generation availability"] = availability

    # Fuel price
    # Load all scenario
     # $/MWh
    
    if settings["Gas price uncertainty"]
        inputs["Variable costs"] = Array(var_cost_input[1:end, 3:end])
        inputs["Gas price scenario probabilities"] = Array(probabilities_input_fuel[probabilities_input_fuel[!,:Uncertainty] .== "Fuel",2:end][1,:])
    else
        inputs["Variable costs"] = Array(var_cost_input[1:end, 2])
        # Deterministic case uses expected fuel price
        inputs["Gas price scenario probabilities"] = [1]
    end
    
    # Keep track of number of scenarios
    inputs["Number of demand scenarios"] = length(inputs["Demand scenario probabilities"])
    inputs["Number of gas price scenarios"] = length(inputs["Gas price scenario probabilities"])
    inputs["Number of weather scenarios"] = length(inputs["Weather scenario probabilities"])
    
    # Storage
    if storage_flag
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
    period_weights_input = CSV.read(string(inputs_path,"/","Representative_period_weights.csv"), DataFrame, header=true)

    n_periods = size(period_weights_input)[1]
    inputs["Number of representative periods"] = n_periods
    T_length = size(demand_input)[1]
    inputs["Number of periods"] = T_length

    # Representative period weights
    t_weights = zeros(T_length, inputs["Number of weather scenarios"])
    if n_periods != T_length
        # Representative periods
        period_length = Integer(T_length/n_periods)
        first_periods = collect(1:period_length:T_length) 
        last_periods = collect(period_length:period_length:T_length)
        inputs["Set of first periods"] = first_periods
        inputs["Set of last periods"] = last_periods

        for scen in 1:inputs["Number of weather scenarios"], i in 1:n_periods
            col = string("Weight-",scen)
            t_weights[first_periods[i]:last_periods[i],scen] .= period_weights_input[!,col][i]/period_length
        end
    else
        # One-to-one mapping
        for scen in 1:inputs["Number of weather scenarios"]
            # col = string("Weight-",scen)
            t_weights[:,scen] = period_weights_input[:,scen+2]
        end
        
    end
    
    inputs["Period weights"] = t_weights
    inputs["Investment costs"] = cost_inv
    inputs["CO2 emission intensities"] = co2_factors
    inputs["Storage flag"] = storage_flag
    inputs["Storage"] = storage_input

    return inputs
end

