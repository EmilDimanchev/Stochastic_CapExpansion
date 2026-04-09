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


####TODO - modify containers to be preallocated sizes to avoid memory over runs from push! and append!

function benders_algorithm(inputs::Dict, settings::Dict, MP::Model, SPs::Array{Model, 3},case_name::String; Eval_SPs = nothing)
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

    conv_tol = settings["Convergence tolerance"]
    capacity_mix = []
    results = Dict{String, Any}()
    output_mp = Dict{String, Any}()
    output_mp_unst = Dict{String, Any}()

    if Eval_SPs !== nothing
        @info("Evaluation SPs provided. Will use these to catalogue performance in addition to base SPs.")
        (S_eval, F_eval, K_eval) = size(Eval_SPs)
        P_eval = inputs["Full Demand scenario probabilities"]
        P_f_eval = inputs["Full Gas price scenario probabilities"]
        P_k_eval = inputs["Full Weather scenario probabilities"]
    end



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
    expected_value = 0.0
    cvar_estimate = 0.0
    min_UB = Inf
    unst_LB = -Inf

    gaps = []

    
    output_mp = run_planning_model(MP, settings)
    alpha_ev = output_mp["Expected alpha"]/settings["Scaling factor cost"]

    if risk_aversion_flag
        u_cvar = output_mp["CVaR term"]/settings["Scaling factor cost"]
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
        sp_all_results = Array{Any}(undef, S, F, K)


        set_capacity_parameters!(SPs, capacity_mix[j])

    
            
        outputs_sp = run_all_subproblems(SPs, inputs, settings, capacity_mix[j])

        sp_all_results = outputs_sp

        sp_obj_per_iter = reshape([outputs_sp[s,f,k]["SP objective"] for s in 1:S, f in 1:F, k in 1:K], (S, F, K))
        
        duals_sp_per_iter =reshape([outputs_sp[s,f,k]["SP dual"][r] for r in 1:R, s in 1:S, f in 1:F, k in 1:K],(R, S, F, K))
        coeffs = reshape([outputs_sp[s,f,k]["coeff"] for s in 1:S, f in 1:F, k in 1:K],(S, F, K))
        expected_value = sum(P[s]*P_f[f]*P_k[k]*sp_obj_per_iter[s,f,k] for s in 1:S, f in 1:F, k in 1:K)/settings["Scaling factor cost"] 
        cvar_estimate = 0.0
    
        if all(size(SPs) .>= (1,1,1))
            cvar_estimate = compute_cvar(sp_obj_per_iter, P, P_f, P_k, VaR_Percent)/settings["Scaling factor cost"]
        end
   
        println("CVaR estimate: ", cvar_estimate)
        println("CVaR Master problem term: ", output_mp["CVaR term"]/settings["Scaling factor cost"])
        if risk_aversion_flag
            UB = risk_aversion_weight*expected_value + (1-risk_aversion_weight)*cvar_estimate + output_mp["Inv_cost"]/settings["Scaling factor cost"]
        else    
            UB = expected_value + output_mp["Inv_cost"]/settings["Scaling factor cost"]
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
            sp_all_results = run_all_subproblems(SPs, inputs, settings, output_mp["Capacity"])

            results["MP"] = output_mp
            results["SPs"] = sp_all_results
            results["CVaR"] = cvar_estimate*settings["Scaling factor cost"]
            results["Expected Value"] = expected_value*settings["Scaling factor cost"]
            results["Capacity per iteration"] = [round.(row; digits=2) for row in capacity_mix]
            #results["Upper bounds"] = upper_bounds
            #results["Lower bounds"] = lower_bounds
            results["LB_hist"] = LB_hist
            results["UB_hist"] = UB_hist
            results["Gaps"] = gaps

            # Evaluate MP sol on Eval SPs if needed
            if Eval_SPs !== nothing
                @info("Evaluating current MP solution on evaluation SPs to track performance on these scenarios.")
                set_capacity_parameters!(Eval_SPs, output_mp["Capacity"])
                outputs_sp_eval = run_all_subproblems(Eval_SPs, inputs, settings, output_mp["Capacity"])
                results["SPs_eval"] = outputs_sp_eval
                results["CVaR"] = compute_cvar(reshape([outputs_sp_eval[s,f,k]["SP objective"] for s in 1:S_eval, f in 1:F_eval, k in 1:K_eval], (S_eval, F_eval, K_eval)), P_eval, P_f_eval, P_k_eval, VaR_Percent)
                results["Expected Value"] = sum(P_eval[s]*P_f_eval[f]*P_k_eval[k]*outputs_sp_eval[s,f,k]["SP objective"] for s in 1:S_eval, f in 1:F_eval, k in 1:K_eval)
            
            end

            break
        else
            add_optimality_cuts!(MP, sp_obj_per_iter, duals_sp_per_iter, capacity_mix[j], coeffs, inputs,settings, j, case_name)
            @info("Running investment problem")
            output_mp_unst = run_planning_model(MP, settings)
            alpha_ev = output_mp_unst["Expected alpha"]/settings["Scaling factor cost"]
            if risk_aversion_flag
                u_cvar = output_mp_unst["CVaR term"]/settings["Scaling factor cost"]
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
                if risk_aversion_flag
                    u_cvar_stb = output_mp["CVaR term"]/settings["Scaling factor cost"]
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
            
        end
    end
    elapsed = time() - algorithm_start_time
    @info("Total elapsed time for algorithm: $(elapsed) seconds")

    # Report results
    @info("Final capacity mix:" * string(results["Capacity per iteration"]))
    return results

end

function compute_cvar(SP_obj, P, P_f, P_k, VaR_Percent)
    alpha = VaR_Percent
    if !(0.0 < alpha <= 1.0)
        throw(ArgumentError("VaR_Percent must be in (0, 1]. Got $(alpha)."))
    end

    # Flatten SP_obj and associated probabilities
    sp_obj_flat = vec(SP_obj)
    prob_flat = vec([
        P[s] * P_f[f] * P_k[k]
        for s in 1:size(SP_obj, 1), f in 1:size(SP_obj, 2), k in 1:size(SP_obj, 3)
    ])

    if length(sp_obj_flat) != length(prob_flat)
        throw(ArgumentError("SP_obj and probability vectors have inconsistent sizes."))
    end

    total_prob = sum(prob_flat)
    if total_prob <= 0.0
        throw(ArgumentError("Total probability mass must be positive."))
    end

    # Normalize probabilities to avoid floating-point drift from 1.0.
    prob_flat = prob_flat ./ total_prob

    # Sort outcomes descending so we integrate directly over the worst alpha tail.
    sorted_indices = sortperm(sp_obj_flat; rev = true)
    sp_obj_sorted = sp_obj_flat[sorted_indices]
    prob_sorted = prob_flat[sorted_indices]

    remaining_tail_mass = alpha
    tail_weighted_sum = 0.0

    for i in eachindex(sp_obj_sorted)
        if remaining_tail_mass <= 1e-14
            break
        end

        weight = min(prob_sorted[i], remaining_tail_mass)
        tail_weighted_sum += sp_obj_sorted[i] * weight
        remaining_tail_mass -= weight
    end

    return tail_weighted_sum / alpha
end

function mga_benders(inputs::Dict, settings::Dict, MP::Model, SPs::Array{Model, 3}, budgets::Dict, case_name::String; Eval_SPs::Union{Array{Model,3}, Nothing} = nothing)
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

    if Eval_SPs !== nothing
        @info("Evaluation SPs provided. Will use these to catalogue performance in addition to base SPs.")
        (S_eval, F_eval, K_eval) = size(Eval_SPs)
        P_eval = inputs["Full Demand scenario probabilities"]
        P_f_eval = inputs["Full Gas price scenario probabilities"]
        P_k_eval = inputs["Full Weather scenario probabilities"]
    end


    conv_tol = settings["Convergence tolerance"]
    capacity_mix = []
    results = Dict{String, Any}()
    output_mp = Dict{String, Any}()
    output_mp_unst = Dict{String, Any}()


    # Initializing new bound settings
    LB_hist = []
    UB_hist = []


    budget_types = collect(keys(budgets))
    UBs = [Inf for _ in budget_types]
    min_UBs = [Inf for _ in budget_types]
    gaps = [Inf for _ in budget_types]
    LBs = [budgets[type] for type in budget_types]
    
    P = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]
    VaR_Percent = settings["Value-at-Risk percent"] 
    risk_aversion_weight = settings["Risk aversion weight"]
    risk_aversion_flag = settings["Risk aversion flag"]
    expected_value = 0.0
    cvar_estimate = 0.0
    gap_hist = []
    

    

    
    output_mp = run_planning_model(MP, settings)

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
        # Actual sps
        @info("Running economic dispatch sub-problems")
        sp_obj_per_iter = zeros(S, F, K)
        duals_sp_per_iter = zeros(R, S, F, K)
        sp_all_results = Array{Any}(undef, S, F, K)


        set_capacity_parameters!(SPs, capacity_mix[j])

    
            
        outputs_sp = run_all_subproblems(SPs, inputs, settings, capacity_mix[j])

        sp_all_results = outputs_sp

        sp_obj_per_iter = reshape([outputs_sp[s,f,k]["SP objective"] for s in 1:S, f in 1:F, k in 1:K], (S, F, K))
        
        duals_sp_per_iter =reshape([outputs_sp[s,f,k]["SP dual"][r] for r in 1:R, s in 1:S, f in 1:F, k in 1:K],(R, S, F, K))
        coeffs = reshape([outputs_sp[s,f,k]["coeff"] for s in 1:S, f in 1:F, k in 1:K],(S, F, K))

        


        #=============

        Compute convergence metrics and update bounds

        ============#
        # Make UB cost estimates based on current SP solutions
        expected_value = sum(P[s]*P_f[f]*P_k[k]*sp_obj_per_iter[s,f,k] for s in 1:S, f in 1:F, k in 1:K)/settings["Scaling factor cost"] 
        cvar_estimate = 0.0
        if isnothing(Eval_SPs) || all(size(SPs) .>= (1,1,1))
            cvar_estimate = compute_cvar(sp_obj_per_iter, P, P_f, P_k, VaR_Percent)/settings["Scaling factor cost"]
        end

        # Build actual UBs based on risk aversion settings and budget types ---- supports multiple budget types at once, but will track the minimum UB across iterations for each type separately
        for (i, budget_type) in enumerate(budget_types)
            if budget_type == "System_Expected"
                UBs[i] = expected_value + output_mp["Inv_cost"]/settings["Scaling factor cost"]
            elseif budget_type == "System_Weighted_CVaR"
                UBs[i] = risk_aversion_weight*expected_value + (1-risk_aversion_weight)*cvar_estimate + output_mp["Inv_cost"]/settings["Scaling factor cost"]
            else
                error("Budget type $(budget_type) not recognized for UB calculation.")
            end

            if UBs[i] < min_UBs[i]
                min_UBs[i] = UBs[i]
            end
        end
        push!(UB_hist, min_UBs)

        
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # MARK: Convergence check
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        for (i, type) in enumerate(budget_types)
            gaps[i] = (min_UBs[i] - LBs[i])/abs(LBs[i])
            @info("Budget type: $(type), Gap: $(gaps[i] * 100)%, LB: $(LBs[i]), UB: $(min_UBs[i])")
        end
        push!(gap_hist, deepcopy(gaps))
        #if maximum(abs.(gaps)) <= conv_tol
        if isapprox(maximum(gaps), 0; atol = settings["Budget convergence tolerance"]) || all(gaps .<= conv_tol)    
            @info("Convergence achieved with maximum percentage gap of $((maximum(gaps)) * 100)%")

            output_mp = run_planning_model(MP, settings)
            set_capacity_parameters!(SPs, output_mp["Capacity"])
            sp_all_results = run_all_subproblems(SPs, inputs, settings, output_mp["Capacity"])

            results["MP"] = output_mp
            results["SPs"] = sp_all_results
            results["Expected Value"] = expected_value*settings["Scaling factor cost"]
            results["Capacity per iteration"] = [round.(row; digits=2) for row in capacity_mix]
            results["CVaR"] = cvar_estimate*settings["Scaling factor cost"]
            #results["Upper bounds"] = upper_bounds
            #results["Lower bounds"] = lower_bounds
            results["UB_hist"] = UB_hist
            results["Gaps"] = gap_hist
            println("Gaps at convergence: ", gap_hist[end])

            # Evaluate MP sol on Eval SPs if needed
            if Eval_SPs !== nothing
                @info("Evaluating current MP solution on evaluation SPs to track performance on these scenarios.")
                set_capacity_parameters!(Eval_SPs, output_mp["Capacity"])
                outputs_sp_eval = run_all_subproblems(Eval_SPs, inputs, settings, output_mp["Capacity"])
                results["SPs_eval"] = outputs_sp_eval
                results["CVaR"] = compute_cvar(reshape([outputs_sp_eval[s,f,k]["SP objective"] for s in 1:S_eval, f in 1:F_eval, k in 1:K_eval], (S_eval, F_eval, K_eval)), P_eval, P_f_eval, P_k_eval, VaR_Percent)
                results["Expected Value"] = sum(P_eval[s]*P_f_eval[f]*P_k_eval[k]*outputs_sp_eval[s,f,k]["SP objective"] for s in 1:S_eval, f in 1:F_eval, k in 1:K_eval)
            end
            break
        else
            add_optimality_cuts!(MP, sp_obj_per_iter, duals_sp_per_iter, capacity_mix[j], coeffs, inputs,settings, j, case_name)
            @info("Running investment problem")
            output_mp = run_planning_model(MP, settings)

            # Update MP solution outputs
            push!(capacity_mix, output_mp["Capacity"])
            @info("New capacity mix: $(round.(output_mp["Capacity"]; digits=2))")
        end
    end
    elapsed = time() - algorithm_start_time
    @info("Total elapsed time for algorithm: $(elapsed) seconds")

    # Report results
    @info("Final capacity mix:" * string(results["Capacity per iteration"]))
    return results

end