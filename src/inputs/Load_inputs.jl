function load_input_data(inputs_path, settings)

    # Collect data here
    inputs = Dict{String, Any}()

    lines_input = CSV.read(string(inputs_path,"/","Network.csv"), DataFrame, header=true)
    inputs["Number of lines"] = size(lines_input)[1]
    # Collect a list of column names matching z1 z2 etc
    zone_columns = filter(x -> occursin(r"^z\d+$", x), names(lines_input))
    inputs["Number of zones"] = length(zone_columns)
    inputs["Line capacities"] = Array(lines_input[:, "Existing_capacity"])
    inputs["Max expansion"] = Array(lines_input[:, "Max_expansion"])
    inputs["CAPEX per MW"] = Array(lines_input[:, "CAPEX_per_MW"])
    inputs["Zone map"] = Array(lines_input[:, zone_columns])

    fuel_response_input = CSV.read(string(inputs_path,"/","Fuel_response.csv"), DataFrame, header=true)
    inputs["Fuel response"] = Array(fuel_response_input[1, 2:end])

    # Technologies
    resources_input = CSV.read(string(inputs_path,"/","Resources.csv"), DataFrame, header=true)

    gen_resources = resources_input[1:end,1]
    inputs["Generation resources"] = gen_resources
    x_max = resources_input[1:end,"Capacity_UB"]
    resource_names = copy(gen_resources)

    # Generation zone map: G × Z binary matrix, entry [g,z] = 1 if resource g is in zone z.
    # Uses the "Zone" column in Resources.csv when present; falls back to all resources in zone 1.
    G = size(resources_input)[1]
    Z = inputs["Number of zones"]
    if "Zone" in names(resources_input)
        gen_zone_map = zeros(Int, G, Z)
        for g in 1:G
            z = resources_input[g, "Zone"]
            gen_zone_map[g, z] = 1
        end
    else
        # Single-zone: all resources belong to zone 1
        gen_zone_map = zeros(Int, G, Z)
        gen_zone_map[:, 1] .= 1
    end

    storage_input = CSV.read(string(inputs_path,"/","Resources_storage.csv"), DataFrame, header=true)
    storage_flag = size(storage_input)[1] > 0

    # Storage zone map: O × Z, then combined with gen into full R × Z resource zone map
    O = size(storage_input)[1]
    if O > 0 && "Zone" in names(storage_input)
        stor_zone_map = zeros(Int, O, Z)
        for o in 1:O
            z = storage_input[o, "Zone"]
            stor_zone_map[o, z] = 1
        end
    else
        # No storage, or single-zone storage: assign all to zone 1
        stor_zone_map = zeros(Int, O, Z)
        if O > 0
            stor_zone_map[:, 1] .= 1
        end
    end
    # Full resource zone map: R × Z  (generation rows first, then storage rows)
    inputs["Resource zone map"] = vcat(gen_zone_map, stor_zone_map)

    var_cost_input = CSV.read(string(inputs_path,"/","Variable_cost.csv"), DataFrame, header=true)
    
    inputs["Number of generation resources"] = size(resources_input)[1]
    inputs["Number of storage resources"] = size(storage_input)[1]
    inputs["Number of resources"] = inputs["Number of generation resources"] + inputs["Number of storage resources"]


    # Stochasticities
    probabilities_input_demand = CSV.read(string(inputs_path,"/","Scenario_probabilities_demand.csv"), DataFrame, header=true)
    probabilities_input_fuel = CSV.read(string(inputs_path,"/","Scenario_probabilities_fuel.csv"), DataFrame, header=true)
    probabilities_input_weather = CSV.read(string(inputs_path,"/","Scenario_probabilities_weather.csv"), DataFrame, header=true)

    # Demand
    demand_scen_files = sort(
        filter(f -> occursin(r"^Demand_\d+\.csv$", basename(f)),
               readdir(inputs_path, join=true)),
        by = f -> parse(Int, match(r"Demand_(\d+)\.csv", basename(f))[1])
    )
    if isempty(demand_scen_files)
        error("No Demand_N.csv files found in $inputs_path")
    end
    demand_scen_data = map(demand_scen_files) do f
        df = CSV.read(f, DataFrame, header=true)

        zone_cols = filter(c -> occursin(r"^Demand(-z\d+)?$", c), names(df))
        Matrix{Float64}(df[:, zone_cols])   # T × Z matrix for this scenario
    end
    
    
    time_index = Array(CSV.read(demand_scen_files[1], DataFrame, header=true)[:, 1])
    
    inputs["Time index"] = time_index
    inputs["Number of periods"] = size(time_index)[1]

    demand_add_input = CSV.read(string(inputs_path,"/","Demand_adders.csv"), DataFrame, header=true)

    inputs["Demand"] = round.(sum(demand_scen_data) ./ length(demand_scen_data), digits=1) #### Why is this the average demand profile?

    

    # Demand shift uncertainty
    if settings["Demand uncertainty"]
        inputs["Demand scenario probabilities"] = Array(probabilities_input_demand[probabilities_input_demand[!,:Uncertainty] .== "Demand",2:end][1,:])
        demand_add = Array(demand_add_input[1,2:end])
    else
        inputs["Demand scenario probabilities"] = [1]
        demand_add = [0]
    end
    inputs["Demand adders"] = demand_add
    
    # Weather uncertainty
    inputs["Demand shift weather"] = zeros(inputs["Number of periods"], Z, 1)
    if settings["Weather uncertainty"]
        # demand_scen_data[k] is T × Z; shift = per-zone deviation from central (average) demand
        # Store as T × Z × K so downstream can index [t, z, k] or [t, k] for single-zone
        K_weather = length(demand_scen_data)
        demand_shift_matrix = zeros(inputs["Number of periods"], Z, K_weather)
        for k in 1:K_weather
            demand_shift_matrix[:, :, k] = round.(demand_scen_data[k] .- inputs["Demand"], digits=1)
        end
        inputs["Demand shift weather"] = demand_shift_matrix

        inputs["Weather scenario probabilities"] = Array(probabilities_input_weather[probabilities_input_weather[!,:Uncertainty] .== "Weather",2:end][1,:])
        
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

        inputs["Demand shift weather"] = zeros(inputs["Number of periods"], Z, 1)
        inputs["Weather scenario probabilities"] = [1]
    end
    
    inputs["Generation availability"] = availability

    # Fuel price
    # Load all scenario
     # $/MWh
    
    if settings["Gas price uncertainty"]
        inputs["Variable costs"] = Array(var_cost_input[1:end, 4:end])
        inputs["Gas price scenario probabilities"] = Array(probabilities_input_fuel[probabilities_input_fuel[!,:Uncertainty] .== "Fuel",2:end][1,:])
    else
        inputs["Variable costs"] = Array(var_cost_input[1:end, 2])
        # Deterministic case uses expected fuel price
        inputs["Gas price scenario probabilities"] = [1]
    end

    if settings["CRM flag"]
        crm_system = CSV.read(string(inputs_path,"/","CRM_System.csv"), DataFrame, header=true)
        crm_resources = CSV.read(string(inputs_path,"/","CRM_Resources.csv"), DataFrame, header=true)

        inputs["CRM Margin"] = Array(crm_system[:, "Margin"])
        inputs["CRM Price Cap"] = Array(crm_system[:, "PriceCap"])
        inputs["CRM Derating Factors"] = Dict(crm_resources[:, "Resource"] .=> crm_resources[:, "Derating_Factor"])
        inputs["Availability CRM"] = mean(inputs["Generation availability"], dims=3)  # Average availability across weather scenarios for CRM derating
        inputs["Demand CRM"] = inputs["Demand"]  # Average demand across scenarios for CRM requirements
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
    
    

    # Representative period weights
    t_weights = zeros(inputs["Number of periods"], inputs["Number of weather scenarios"])
    if n_periods != inputs["Number of periods"]
        # Representative periods
        period_length = Integer(inputs["Number of periods"]/n_periods)
        first_periods = collect(1:period_length:inputs["Number of periods"]) 
        last_periods = collect(period_length:period_length:inputs["Number of periods"])
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

