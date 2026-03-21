# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stochastic equilibrium model algorithm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


export run_benders_algorithm, benders_algorithm

function run_benders_algorithm(inputs::Dict, settings::Dict)

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
            
        @info("Iteration: $(j), Maximum percentage gap: $(maximum(abs.(gaps)) * 100)%")


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
    @info("Final capacity mix: $(results["Capacity per iteration"][end])")
    return results
end


function benders_algorithm(inputs::Dict, settings::Dict, MP::Model, SPs::Array{Model, 3})
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
    MP_obj = []
    conv_tol = settings["Convergence tolerance"]
    lower_bounds = []
    upper_bounds = []
    capacity_mix = []
    cvar = []
    alphas = []
    inv_cost = []
    results = Dict{String, Any}()
    #capacity_mix_initial = ones(R).*5e3
    #push!(capacity_mix, capacity_mix_initial)
    push!(lower_bounds, zeros(S,F,K).*Inf)
    push!(cvar, 0)
    output_mp = Dict{String, Any}()
    output_mp_unst = Dict{String, Any}()


    # Initializing new bound settings
    LB_hist = []
    UB_hist = []
    
    UB = Inf
    P = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]
    VaR_Percent = settings["Value-at-Risk percent"] 
    risk_aversion_weight = settings["Risk aversion weight"]
    risk_aversion_flag = settings["Risk aversion flag"]
    expected_value_hist = []
    cvar_hist = []
    alpha_ev_hist = []
    u_cvar_hist = []
    coeffs_hist = []
    min_UB = Inf
    unst_LB = -Inf

    gaps = []

    
    output_mp = run_planning_model(MP, settings)
    alpha_ev = output_mp["Expected alpha"]/settings["Scaling factor cost"]
    push!(alpha_ev_hist, alpha_ev)
    push!(inv_cost, output_mp["Inv_cost"])

    if risk_aversion_flag
        u_cvar = output_mp["CVaR term"]/settings["Scaling factor cost"]
        push!(u_cvar_hist, u_cvar)
        risk_adjusted_LB = risk_aversion_weight*alpha_ev + (1-risk_aversion_weight)*u_cvar
        LB = risk_adjusted_LB + output_mp["Inv_cost"]/settings["Scaling factor cost"]
    else
        LB = alpha_ev + output_mp["Inv_cost"]/settings["Scaling factor cost"]
    end
    capacity_mix_initial = output_mp["Capacity"]
    push!(capacity_mix, capacity_mix_initial)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Algorithm
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    algorithm_start_time = time()

    for j in 1:J_max
        println()
        println(string("*** Iteration: ", j))

        
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


        set_capacity_parameters!(SPs, capacity_mix[j])

    
            
        outputs_sp = run_all_subproblems(SPs, inputs, settings)

        sp_all_results = outputs_sp

        sp_obj_per_iter = reshape([outputs_sp[s,f,k]["SP objective"] for s in 1:S, f in 1:F, k in 1:K], (S, F, K))
        
        duals_sp_per_iter =reshape([outputs_sp[s,f,k]["SP dual"][r] for r in 1:R, s in 1:S, f in 1:F, k in 1:K],(R, S, F, K))
        coeffs = reshape([outputs_sp[s,f,k]["coeff"] for s in 1:S, f in 1:F, k in 1:K],(S, F, K))
        push!(coeffs_hist, coeffs)
        
        # Update SP solution outputs 
        push!(SP_obj, sp_obj_per_iter)
        push!(SP_duals, duals_sp_per_iter)

        push!(expected_value_hist, sum(P[s]*P_f[f]*P_k[k]*sp_obj_per_iter[s,f,k] for s in 1:S, f in 1:F, k in 1:K)/settings["Scaling factor cost"]) 
        cvar_estimate = compute_cvar(SP_obj[j], P, P_f, P_k, VaR_Percent)/settings["Scaling factor cost"]
        if risk_aversion_flag
            push!(cvar_hist, cvar_estimate)
            UB = risk_aversion_weight*expected_value_hist[end] + (1-risk_aversion_weight)*cvar_estimate + output_mp["Inv_cost"]/settings["Scaling factor cost"]
        else    
            UB = expected_value_hist[end] + output_mp["Inv_cost"]/settings["Scaling factor cost"]
        end
        if UB < min_UB
            min_UB = UB
        end
        push!(UB_hist, min_UB)
        
        
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # MARK: Convergence check
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # UB
        #ub_estimate = SP_obj[j]./settings["Scaling factor cost"]
        #push!(upper_bounds, ub_estimate)
        
        #gaps = (upper_bounds[j].-lower_bounds[j])./lower_bounds[j]  
        
        #if any(1e-6 .< gaps .< 0)
        #    @warn("Negative gap detected; check model formulation.")
        #end

            
        #@info("Maximum percentage gap: $(maximum(abs.(gaps)) * 100)%")

        gap = (min_UB - LB)/abs(LB)
        push!(gaps, gap)
        @info("Gap: $(gap * 100)%")
        @info("LB: $(LB), UB: $(min_UB)")
        if gap < 0 && abs(gap) > 1e-6
            @warn("Negative gap detected; check model formulation.")
        end


        #if maximum(abs.(gaps)) <= conv_tol
        if gap <= conv_tol    
            @info("Convergence achieved with maximum percentage gap of $((gap) * 100)%")

            output_mp = run_planning_model(MP, settings)
            set_capacity_parameters!(SPs, output_mp["Capacity"])
            sp_all_results = run_all_subproblems(SPs, inputs, settings)

            results["MP"] = output_mp
            results["SPs"] = sp_all_results
            results["CVaR"] = cvar_estimate*settings["Scaling factor cost"]
            results["Expected Value"] = expected_value_hist[end]*settings["Scaling factor cost"]
            results["Capacity per iteration"] = [round.(row; digits=2) for row in capacity_mix]
            #results["Upper bounds"] = upper_bounds
            #results["Lower bounds"] = lower_bounds
            results["LB_hist"] = LB_hist
            results["UB_hist"] = UB_hist

            break
        else
            add_optimality_cuts!(MP, SP_obj[j], SP_duals[j], capacity_mix[j], coeffs, inputs,settings, j)
            @info("Running investment problem")
            output_mp_unst = run_planning_model(MP, settings)
            alpha_ev = output_mp_unst["Expected alpha"]/settings["Scaling factor cost"]
            #push!(alpha_ev_hist, alpha_ev)
            #push!(inv_cost, output_mp["Inv_cost"])
            if risk_aversion_flag
                u_cvar = output_mp_unst["CVaR term"]/settings["Scaling factor cost"]
                #push!(u_cvar_hist, u_cvar)
                risk_adjusted_LB = risk_aversion_weight*alpha_ev + (1-risk_aversion_weight)*u_cvar
                LB_unst = risk_adjusted_LB + output_mp_unst["Inv_cost"]/settings["Scaling factor cost"]
            else
                LB_unst  = alpha_ev + output_mp_unst["Inv_cost"]/settings["Scaling factor cost"]
            end
            LB = max(LB, LB_unst)
            push!(LB_hist, LB)
            
            if settings["Regularization flag"]
                @info("Applying regularization to master problem")
                output_mp = regularization(MP, min_UB, LB, settings)
                alpha_ev_stb = output_mp["Expected alpha"]/settings["Scaling factor cost"]
                #push!(alpha_ev_hist, alpha_ev)
                #push!(inv_cost, output_mp["Inv_cost"])
                if risk_aversion_flag
                    u_cvar_stb = output_mp["CVaR term"]/settings["Scaling factor cost"]
                    #push!(u_cvar_hist, u_cvar)
                    risk_adjusted_LB_stb = risk_aversion_weight*alpha_ev_stb + (1-risk_aversion_weight)*u_cvar_stb
                    LB = risk_adjusted_LB_stb + output_mp["Inv_cost"]/settings["Scaling factor cost"]
                else
                    LB = alpha_ev_stb + output_mp["Inv_cost"]/settings["Scaling factor cost"]
                end 
            else
                output_mp = output_mp_unst
            end

            
        

            # Update MP solution outputs
            push!(capacity_mix, output_mp["Capacity"])
            @info("New capacity mix: $(round.(output_mp["Capacity"]; digits=2))")
            push!(MP_obj, output_mp["Planning objective"])
            push!(alphas, output_mp["Alpha"])
            push!(inv_cost, output_mp["Inv_cost"])
            
            
            # LB
            push!(lower_bounds, output_mp["Alpha"]./settings["Scaling factor cost"])

            alpha_ev = output_mp["Expected alpha"]/settings["Scaling factor cost"]
            push!(alpha_ev_hist, alpha_ev)
            push!(inv_cost, output_mp["Inv_cost"])
            if risk_aversion_flag
                u_cvar = output_mp["CVaR term"]/settings["Scaling factor cost"]
                push!(u_cvar_hist, u_cvar)
                risk_adjusted_LB = risk_aversion_weight*alpha_ev + (1-risk_aversion_weight)*u_cvar
                LB = risk_adjusted_LB + output_mp["Inv_cost"]/settings["Scaling factor cost"]
            else
                LB = alpha_ev + output_mp["Inv_cost"]/settings["Scaling factor cost"]
            end
            
        end
    end
    elapsed = time() - algorithm_start_time
    @info("Total elapsed time for algorithm: $(elapsed) seconds")

    # Report results
    @info("Final capacity mix:" * string(results["Capacity per iteration"]))
    return results, gaps

end

function compute_cvar(SP_obj, P, P_f, P_k, VaR_Percent)
    # Flatten SP_obj and associated probabilities
    sp_obj_flat = vec(SP_obj)
    prob_flat = vec([P[s]*P_f[f]*P_k[k] for s in 1:size(SP_obj,1), f in 1:size(SP_obj,2), k in 1:size(SP_obj,3)])
    
    # Sort SP outcomes and probabilities
    sorted_indices = sortperm(sp_obj_flat)
    sp_obj_sorted = sp_obj_flat[sorted_indices]
    prob_sorted = prob_flat[sorted_indices]
    
    # Compute cumulative probabilities
    cum_prob = cumsum(prob_sorted)
    
    # Identify VaR threshold
    var_threshold_index = findfirst(x -> x >= (1-VaR_Percent)+.005, cum_prob) # Adding a small tolerance to handle floating-point issues

    # Compute CVaR as weighted average of tail outcomes
    cvar = sum(sp_obj_sorted[i]*prob_sorted[i] for i in var_threshold_index:length(sp_obj_sorted))/(VaR_Percent)
    
    return cvar
end