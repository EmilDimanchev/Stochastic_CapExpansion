# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stochastic equilibrium model algorithm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using Dates

# Utilities
include("./utilities/write_results.jl")
include("./utilities/load_inputs.jl")

# Load models
include("./subproblems_economic_dispatch.jl")
include("./master_planning.jl")

# Set up experiment
include("./create_case.jl")

# Load data
inputs = load_input_data(inputs_path, settings)

# Iterations
J_max = 150

# Initialize
settings["Search equilibria"] = false
S = inputs["Number of demand scenarios"]
F = inputs["Number of gas price scenarios"]
K = inputs["Number of weather scenarios"]
G = inputs["Number of generation resources"]    
O =  inputs["Number of storage resources"]
R = G + O
T = inputs["Number of periods"]
SP_obj = []
SP_duals = []
power_price = []
nse_dual_cost = []
consumer_surplus = []
MP_obj = []
conv_tol = settings["Convergence tolerance"]
lower_bounds = []
upper_bounds = []
capacity_mix = []
cvar = []
alphas = []
results = Dict{String, Any}()
capacity_mix_initial = ones(R).*5e3
push!(capacity_mix, capacity_mix_initial)
push!(lower_bounds, zeros(S,F,K).*Inf)
push!(cvar, 0)
one_last_iter_flag = false

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Algorithm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

algorithm_start_time = time()

for j in 1:J_max
    
    @info(string("*** Iteration: ", j))
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MARK: Master
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    @info("Running investment problem")
    if j >= 2
      
      
        output_mp = run_planning_model(inputs, settings, SP_obj, SP_duals, capacity_mix)

        # Update MP solution outputs
        push!(capacity_mix, output_mp["Capacity"])
        println("New capacity mix: ", round.(output_mp["Capacity"]; digits=2))
        push!(MP_obj, output_mp["Planning objective"])
        push!(alphas, output_mp["Alpha"])
        
        # LB
        push!(lower_bounds, output_mp["Alpha"]./settings["Scaling factor cost"])
        
    end
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MARK: Subproblems
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    @info("Running economic dispatch sub-problems")
    sp_obj_per_iter = zeros(S, F, K)
    duals_sp_per_iter = zeros(R, S, F, K)
    power_price_per_iter = zeros(T, S, F, K)
    nse_dual_cost_per_iter = zeros(S, F, K)
    consumer_surplus_per_iter = zeros(S, F, K)
    sp_all_results = Array{Any}(undef, S, F, K)

    for s in 1:S, f in 1:F, k in 1:K

        output_sp = run_economic_dispatch(inputs, settings, capacity_mix[j], inputs["Demand shift weather"][:,k], inputs["Variable costs"][:,f], inputs["Generation availability"][:,:,k], inputs["Period weights"][:,k], inputs["Demand adders"][s])

        sp_all_results[s,f,k] = output_sp

        sp_obj_per_iter[s,f,k] = output_sp["SP objective"]
        
        duals_sp_per_iter[:,s,f,k] = output_sp["SP dual"]
        
    end
    
    # Update SP solution outputs 
    push!(SP_obj, sp_obj_per_iter)
    push!(SP_duals, duals_sp_per_iter)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MARK: Convergence check
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # UB
    ub_estimate = SP_obj[j]./settings["Scaling factor cost"]
    push!(upper_bounds, ub_estimate)
    
    gaps = (upper_bounds[j].-lower_bounds[j])./lower_bounds[j]  
    
    if any(1e-6 .< gaps .< 0)
        @warn("Negative gap detected; check model formulation.")
    end

    if maximum(abs.(gaps)) <= conv_tol
        
        @info("Convergence achieved with maximum percentage gap of $(maximum(abs.(gaps)) * 100)%")

        results["MP"] = output_mp
        results["SP"] = sp_all_results
        results["Capacity per iteration"] = [round.(row; digits=2) for row in capacity_mix]
        results["Upper bounds"] = upper_bounds
        results["Lower bounds"] = lower_bounds

        break
    end
end
elapsed = time() - algorithm_start_time
@info("Total elapsed time for algorithm: $(elapsed) seconds")

# Report results
total_ub = [sum(results["Upper bounds"][i]) for i in 1:length(upper_bounds)]
total_lb = [sum(results["Lower bounds"][i]) for i in 1:length(lower_bounds)]
total_gap = total_ub .- total_lb
perc_gap = total_gap ./ total_lb
results["Capacity per iteration"]