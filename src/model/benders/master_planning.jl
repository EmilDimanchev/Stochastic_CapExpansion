export run_planning_model, build_planning_model, add_optimality_cuts!

#=====================================

Base function (Out of date, but keeping for reference)

=====================================#

function run_planning_model(inputs, settings, SP_obj_list, SP_dual_list, x_prev)
    # Model
    MP = Model(Gurobi.Optimizer)
    set_silent(MP)
    set_optimizer_attribute(MP, "OptimalityTol", 1e-5)
    set_optimizer_attribute(MP, "FeasibilityTol", 1e-5)
    set_optimizer_attribute(MP, "Crossover", 0) # use automatic crossover
    
    # Sets and flags
    risk_aversion_flag = settings["Risk aversion flag"]
    risk_aversion_weight = settings["Risk aversion weight"]

    cost_inv = inputs["Investment costs"]
    P = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]

    scaling_factor_cost = settings["Scaling factor cost"]
    cost_inv = cost_inv./scaling_factor_cost
    capex_line = inputs["CAPEX per MW"]./scaling_factor_cost

    # ~~~
    # Build model
    # ~~~

    # Sets
    J = length(SP_obj_list)
    S = size(P)[1] # number of demand scenarios
    F = size(P_f)[1] # number of gas price scenarios
    K = size(P_k)[1] # number of weather scenarios
    G = inputs["Number of generation resources"]
    O =  inputs["Number of storage resources"]
    L = inputs["Number of lines"]
    R = G + O


    # ~~~
    # Parameters
    # ~~~

    if risk_aversion_flag
        # CVaR parameters
        Ψ = settings["Value-at-Risk percent"] 
        Ω = risk_aversion_weight
    else
        Ω = 1
    end
    
    # set_optimizer_attribute(MP, "OutputFlag", 0)

    # ~~~~
    # Load inputs
    # ~~~

    SP_obj = reduce((a, b) -> cat(a, b; dims=4), SP_obj_list)
    SP_dual = reduce((a, b) -> cat(a, b; dims=5), SP_dual_list)

    
    SP_obj = SP_obj./scaling_factor_cost
    SP_dual = SP_dual./scaling_factor_cost
    x_ub = 1e6

    # ~~~
    # Model formulation
    # ~~~

    # Capacity
    @variable(MP, x[r in 1:R] >= 0) # Capacity, MW
    @constraint(MP, max_capacity[r in 1:R], x[r] <= x_ub)
    @variable(MP, x_line[l in 1:L] >= 0) # Transmission line expansion, MW
    @constraint(MP, max_line_expansion[l in 1:L], x_line[l] <= inputs["Max expansion"][l])

    # Cuts
    @variable(MP, alpha[s in 1:S, f in 1:F, k in 1:K] >= 0)

    if risk_aversion_flag
        # Auxiliary varliables for CVaR
        @variable(MP, u[s in 1:S, f in 1:F, k in 1:K] >= 0) # loss relative to VaR, $/MW
        @variable(MP, ζ) # VaR variable, $/MW
        @constraint(MP, cvar_tail[s in 1:S, f in 1:F, k in 1:K, j in 1:J], MP[:u][s,f,k] >= SP_obj[s,f,k,j] + sum(SP_dual[r,s,f,k,j]*(MP[:x][r]-x_prev[j][r]) for r in 1:R) - MP[:ζ])
    end
        # @constraint(MP, cvar_tail[s in 1:S, f in 1:F, k in 1:K], u[s,f,k] >= sum(t_weights[t]*g[r,t,s,f,k]*cost_var[r,f] for r in 1:G) + sum(t_weights[t]*eNSE_Cost[t,s,f,k] for t in 1:T) + sum(cost_inv[r]*x[r] for r in 1:R) - ζ)
    


   
    @constraint(MP, optimality_cuts[s in 1:S, f in 1:F, k in 1:K, j in 1:J], MP[:alpha][s,f,k] >= SP_obj[s,f,k,j] + sum(SP_dual[r,s,f,k,j]*(MP[:x][r]-x_prev[j][r]) for r in 1:R))


    # ~~~
    # Objective function
    # ~~~ 

    if risk_aversion_flag
        @objective(MP, Min, sum(x[r]*cost_inv[r] for r in 1:R) + Ω*sum(P[s]*P_f[f]*P_k[k]*alpha[s,f,k] for s in 1:S, f in 1:F, k in 1:K) + (1-Ω)*(ζ + 1/Ψ*sum(P[s]*P_f[f]*P_k[k]*u[s,f,k] for s in 1:S, f in 1:F, k in 1:K)))
    else
        @objective(MP, Min, sum(x[r]*cost_inv[r] for r in 1:R) + sum(P[s]*P_f[f]*P_k[k]*alpha[s,f,k] for s in 1:S, f in 1:F, k in 1:K))
    end


    optimize!(MP)

    # Write outputs
    output = Dict{String, Any}()
    output["Planning objective"] = objective_value(MP).*scaling_factor_cost
    output["Capacity"] = value.(x)
    output["Alpha"] = value.(alpha).*scaling_factor_cost
    output["Contract volume"] = 0
    output["Contract price"] = 0
    
    # Record results from the linearization or not
    if risk_aversion_flag
        output["CVaR Loss"] = value.(u).*scaling_factor_cost
        output["VaR"] = value(ζ)*scaling_factor_cost
    else
        output["CVaR Loss"] = zeros(S)
        output["VaR"] = 0
    end

    return output
end

#=====================================

Constructors

=====================================#

function build_planning_model(inputs, settings)
    # Model
    MP = Model()
    set_silent(MP)

    if settings["Solver"] == "HiGHS"
        set_optimizer(MP, HiGHS.Optimizer)
        set_optimizer_attribute(MP, "solver", "choose")
        set_optimizer_attribute(MP, "run_crossover", "off")
    elseif settings["Solver"] == "Gurobi"
        set_optimizer(MP, Gurobi.Optimizer)
        set_optimizer_attribute(MP, "OptimalityTol", 1e-5)
        set_optimizer_attribute(MP, "FeasibilityTol", 1e-5)
        set_optimizer_attribute(MP, "Crossover", 0)
    end

    # Sets and flags
    risk_aversion_flag = settings["Risk aversion flag"]

    cost_inv = inputs["Investment costs"]
    P = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]

    scaling_factor_cost = settings["Scaling factor cost"]
    cost_inv = cost_inv./scaling_factor_cost
    capex_line = inputs["CAPEX per MW"]./scaling_factor_cost

    # ~~~
    # Build model
    # ~~~

    # Sets
    #J = length(SP_obj_list)
    S = size(P)[1] # number of demand scenarios
    F = size(P_f)[1] # number of gas price scenarios
    K = size(P_k)[1] # number of weather scenarios
    G = inputs["Number of generation resources"]
    O =  inputs["Number of storage resources"]
    R = G + O
    Z = inputs["Number of zones"]
    L = inputs["Number of lines"]

    Ω = settings["Risk aversion weight"]
    x_ub = 1e6

    # ~~~
    # Model formulation
    # ~~~

    # Capacity
    @variable(MP, x[r in 1:R] >= 1e-5) # Capacity, MW
    @constraint(MP, max_capacity[r in 1:R], x[r] <= x_ub)

    @variable(MP, x_line[l in 1:L] >= 0) # Transmission line expansion, MW
    @constraint(MP, max_line_expansion[l in 1:L], x_line[l] <= inputs["Max expansion"][l])

    # Cuts
    @variable(MP, alpha[s in 1:S, f in 1:F, k in 1:K] >= 0)
    

    # Cost Expressions
    @expression(MP, inv_cost, sum(x[r]*cost_inv[r] for r in 1:R) + sum(x_line[l]*capex_line[l] for l in 1:L))
    @expression(MP, expected_alpha, sum(P[s]*P_f[f]*P_k[k]*alpha[s,f,k] for s in 1:S, f in 1:F, k in 1:K))
    @expression(MP, exp_sys_cost, expected_alpha + inv_cost)
    if risk_aversion_flag
        add_risk_terms!(MP, inputs, settings)
    end

    @expression(MP, inv_cost_by_zone[z in 1:Z], sum(x[r]*cost_inv[r]*inputs["Resource zone map"][r,z] for r in 1:R))

    # ~~~
    # Objective function
    # ~~~ 

    if risk_aversion_flag
        @expression(MP, eObj, inv_cost + Ω*expected_alpha + (1-Ω)*MP[:cvar_term])
        @objective(MP, Min, eObj)
    else
        @expression(MP, eObj, inv_cost + expected_alpha)
        @objective(MP, Min, eObj)
    end
    return MP

end

function add_risk_terms!(MP::Model, inputs, settings)
    # Sets and flags
    risk_aversion_weight = settings["Risk aversion weight"]

    P = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]

    S = size(P)[1] # number of demand scenarios
    F = size(P_f)[1] # number of gas price scenarios
    K = size(P_k)[1] # number of weather scenarios

    # CVaR parameters
    Ψ = settings["Value-at-Risk percent"] 
    Ω = risk_aversion_weight

    # Auxiliary varliables for CVaR
    @variable(MP, u[s in 1:S, f in 1:F, k in 1:K] >= 0) # loss relative to VaR, $/MW
    @variable(MP, ζ) # VaR variable, $/MW
    
    # CVaR definition constraint
    @constraint(MP, cvar_tail[s in 1:S, f in 1:F, k in 1:K], MP[:u][s,f,k] >= -MP[:ζ])
    
    # Accounting for risk in the objective function
    @expression(MP, cvar_term, ζ + 1/Ψ*sum(P[s]*P_f[f]*P_k[k]*u[s,f,k] for s in 1:S, f in 1:F, k in 1:K))
    @expression(MP, risk_adjusted_sys_cost, MP[:inv_cost] + Ω*MP[:expected_alpha] + (1-Ω)*MP[:cvar_term])

end

function add_optimality_cuts!(MP, SP_obj, cap_dual, line_dual, x_prev, x_prev_line, coeffs, inputs,settings, iteration, case)
    scaling_factor_cost = settings["Scaling factor cost"]
    cvar_tail = settings["Risk aversion flag"]

    # Input Sets
    P = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]

    # Set Sizes
    S = size(P)[1] # number of demand scenarios
    F = size(P_f)[1] # number of gas price scenarios
    K = size(P_k)[1] # number of weather scenarios
    G = inputs["Number of generation resources"]
    O =  inputs["Number of storage resources"]
    R = G + O
    L = inputs["Number of lines"]

    SP_obj = SP_obj./scaling_factor_cost
    cap_dual = cap_dual./scaling_factor_cost
    line_dual = line_dual./scaling_factor_cost
    if cvar_tail
        @constraint(MP, [s in 1:S, f in 1:F, k in 1:K], coeffs[s,f,k]*MP[:u][s,f,k] >= SP_obj[s,f,k] + sum(cap_dual[r,s,f,k]*(MP[:x][r]-x_prev[r]) for r in 1:R) + sum(line_dual[l,s,f,k]*(MP[:x_line][l]-x_prev_line[l]) for l in 1:L) - MP[:ζ], base_name = "cvar_tail_cuts_"*case*"_"*string(iteration)*"_")
    end
    
    @constraint(MP, [s in 1:S, f in 1:F, k in 1:K], coeffs[s,f,k]*MP[:alpha][s,f,k] >= SP_obj[s,f,k] + sum(cap_dual[r,s,f,k]*(MP[:x][r]-x_prev[r]) for r in 1:R) + sum(line_dual[l,s,f,k]*(MP[:x_line][l]-x_prev_line[l]) for l in 1:L), base_name = "optimality_cut_"*case*"_"*string(iteration)*"_")
end

#=====================================

Runner

=====================================#

# Multiple dispatch run_planning_model with and without SP outputs for cuts
function run_planning_model(MP, settings)

    optimize!(MP)
    if termination_status(MP) != MOI.OPTIMAL
        @warn("Model did not solve to optimality. Status: ", termination_status(MP))
    end
    if termination_status(MP) == MOI.INFEASIBLE_OR_UNBOUNDED || termination_status(MP) == MOI.INFEASIBLE
        @warn("Model did not solve to optimality. Status: ", termination_status(MP))
        compute_conflict!(MP)
        list_of_conflicting_constraints = ConstraintRef[];
        for (F, S) in list_of_constraint_types(MP)
            for con in all_constraints(MP, F, S)
                if get_attribute(con, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
                    push!(list_of_conflicting_constraints, con)
                end
            end
        end
        display(list_of_conflicting_constraints)
        error("Model infeasible. See conflicting constraints above.")
    end

    # Write outputs
    output = write_outputs(MP, settings)

    return output

end


#=====================================

Regularization functions

=====================================#

function adjust_gamma(prev_UB, UBnew, LB, gamma)
    Ia = prev_UB - UBnew;
	Ip = prev_UB - LB-gamma*(prev_UB-LB);
	ω = 0.5;
	η₁ = 0.1;
	η₂ = 0.9;
    max_gamma = 0.25
	if Ia>=0 && Ip>0

		r = Ia/Ip
		
		if r <=η₁
			#bad iteration: there wasn't enough improvement, we increase γ because we are less confident in the piecewise approximation in the master model
			gamma = gamma*1.2;
		elseif η₁ <= r <= η₂
			#okish iteration
			#do nothing
		else
			#good iteration, we are feeling confident in the piecewise approximation in the master model and so we reduce γ
			gamma = gamma/1.2;
		end
        gamma = min(gamma, max_gamma)

	else
		#do nothing
	end
    return gamma
end


function level_set_regularization(MP, UB, LB, gamma, settings)

    @constraint(MP, cLevel_set, MP[:eObj] <= LB + 0.5*(UB-LB))
    #@constraint(MP, cLevel_set, MP[:eObj] >= LB + gamma*(UB-LB))
    if settings["Solver"] == "HiGHS"
        set_optimizer_attribute(MP, "solver", "ipm")
    elseif settings["Solver"] == "Gurobi"
        set_optimizer_attribute(MP, "Method", 2)
    end
    @objective(MP, Min, 0.0)

    optimize!(MP)

    if termination_status(MP) != MOI.OPTIMAL
        @warn("Model did not solve to optimality. Status: ", termination_status(MP))
    end
    if termination_status(MP) == MOI.INFEASIBLE
        @warn("Model did not solve to optimality. Status: ", termination_status(MP))
        compute_conflict!(MP)
        list_of_conflicting_constraints = ConstraintRef[];
        for (F, S) in list_of_constraint_types(MP)
            for con in all_constraints(MP, F, S)
                if get_attribute(con, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
                    push!(list_of_conflicting_constraints, con)
                end
            end
        end
        display(list_of_conflicting_constraints)
        error("Model infeasible. See conflicting constraints above.")
    end
    
    # Write outputs
    output = write_outputs(MP, settings)

    delete(MP,MP[:cLevel_set])
    unregister(MP,:cLevel_set)
    @objective(MP,Min, MP[:eObj])
    if settings["Solver"] == "HiGHS"
        set_optimizer_attribute(MP, "solver", "choose")
    elseif settings["Solver"] == "Gurobi"
        set_optimizer_attribute(MP, "Method", -1)
    end

    return output

end

# Quadratic trust region regularization
# per https://www.sciencedirect.com/science/article/pii/S0377221719303617
function qtr_regularization(MP, stab_cent_cap, stab_cent_line, gamma_qtr, settings)
    # Sets
    R = length(MP[:x])
    L = length(MP[:x_line])

    @constraint(MP, cTrust_region, sum((MP[:x][r]-stab_cent_cap[r])^2 for r in 1:R) + sum((MP[:x_line][l]-stab_cent_line[l])^2 for l in 1:L) <= gamma_qtr)

    optimize!(MP)

    if termination_status(MP) != MOI.OPTIMAL
        @warn("Model did not solve to optimality. Status: ", termination_status(MP))
    end
    if termination_status(MP) == MOI.INFEASIBLE
        @warn("Model did not solve to optimality. Status: ", termination_status(MP))
        compute_conflict!(MP)
        list_of_conflicting_constraints = ConstraintRef[];
        for (F, S) in list_of_constraint_types(MP)
            for con in all_constraints(MP, F, S)
                if get_attribute(con, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
                    push!(list_of_conflicting_constraints, con)
                end
            end
        end
        display(list_of_conflicting_constraints)
        error("Model infeasible. See conflicting constraints above.")
    end
    
    # Write outputs
    output = write_outputs(MP, settings; reg=true)

    delete(MP,MP[:cTrust_region])
    unregister(MP,:cTrust_region)
    return output
end


function set_gamma_qtr(gamma_qtr, phi, U_unst, U_prev, x_vector)
    phi_min = 0.1
    kappa = 0.5
    if abs(1- (U_unst/U_prev)) < 0.01
        phi = max(phi_min, kappa*phi)
    end

    gamma_qtr = phi^2 * norm(x_vector, 1)^2

    return gamma_qtr, phi
end




#=====================================

Outputs

=====================================#


function write_outputs(MP, settings; reg=false)
    scaling_factor_cost = settings["Scaling factor cost"]
    output = Dict{String, Any}()
    output["Planning objective"] = value(MP[:eObj])*scaling_factor_cost
    output["Capacity"] = value.(MP[:x])
    output["Line expansion"] = value.(MP[:x_line])
    output["Alpha"] = value.(MP[:alpha]).*scaling_factor_cost
    output["Contract volume"] = 0
    output["Contract price"] = 0
    output["Inv_cost"] = value(MP[:inv_cost])*scaling_factor_cost
    output["Expected alpha"] = value(MP[:expected_alpha])*scaling_factor_cost
    output["Expected System Cost"] = value.(MP[:exp_sys_cost])*scaling_factor_cost
    output["Inv cost by zone"] = value.(MP[:inv_cost_by_zone]).*scaling_factor_cost
    
    # Record results from the linearization or not
    if settings["Risk aversion flag"]
        output["CVaR Loss"] = value.(MP[:u]).*scaling_factor_cost
        output["VaR"] = value.(MP[:ζ]).*scaling_factor_cost
        output["CVaR term"] = value(MP[:cvar_term])*scaling_factor_cost
        output["Risk adjusted system cost"] = value(MP[:risk_adjusted_sys_cost])*scaling_factor_cost
    else
        output["CVaR Loss"] = 0
        output["VaR"] = 0
        output["CVaR term"] = 0
    end
    
    if !reg
        cuts = [con for con in all_constraints(MP, include_variable_in_set_constraints=false) if startswith(string(con), "optimality_cut") || startswith(string(con), "cvar_tail_cuts")]
        cut_names = [name(con) for con in cuts]
        cut_duals = [dual(con) for con in cuts]
        output["Cut Duals"] = Dict(zip(cut_names, cut_duals))
    end

    return output

end
#=====================================

MGA Functions for Benders Implementation

======================================#

function set_objective_bendersMP!(model, objective_type::String, inputs, settings; obj_weight::Float64 = -1.0, set_coeffs::Vector = [])
    if objective_type == "System_Expected"
        unregister(model, :eObj)
        @expression(model, eObj, model[:exp_sys_cost])
        @objective(model, Min, model[:eObj])
        @info("Objective set to minimize expected system cost: ", objective_function(model))
    elseif objective_type == "System_Weighted_CVaR"
        if obj_weight < -0.0001 || obj_weight > 1.001
            error("For Weighted CVaR objective, please provide a valid weight between 0 and 1.")
        end
        # Check if model is initialized with risk terms?
        if !haskey(model, :u) || !haskey(model, :cvar_term)
            @info("Initializing risk terms for Weighted CVaR objective.")
            add_risk_terms!(model, inputs, settings)
            settings["risk aversion flag"] = true
        end
        unregister(model, :eObj)
        @expression(model, eObj, model[:inv_cost] + (obj_weight)*model[:expected_alpha] + (1-obj_weight)*model[:cvar_term])
        @objective(model, Min, model[:eObj])
        @info("Objective set to minimize weighted CVaR. Weight is: ", obj_weight, " Objective function: ", objective_function(model))
    elseif objective_type == "Capacity"
        if set_coeffs == []
            error("For Capacity objective, please provide a vector of coefficients for the capacity expression.")
        end
        x_full = [collect(model[:x][r] for r in 1:length(model[:x])); collect(model[:x_line][l] for l in 1:length(model[:x_line]))]
        @objective(model, Min, set_coeffs' * x_full)
        @info("Objective set to minimize capacity expression: ", objective_function(model))
    else
        error("Unsupported objective type. Please choose from 'Expectation', 'Weighted CVaR', or 'Capacity'.")
    end
end




function add_budget_constraint_bendersMP(model::Model, budget::Float64, budget_type::String, existing_budgets::Dict; risk_aversion::Float64 = -1.0)
    budget_added = Symbol[]
    if budget_type == "Investment"
        @constraint(model, c_budget_inv, model[:inv_cost] <= budget)
        push!(budget_added, :c_budget_inv)
        existing_budgets[budget_type] = budget
    elseif budget_type == "Expected"
        @constraint(model, c_budget_exp, model[:expected_alpha] <= budget)
        push!(budget_added, :c_budget_exp)
        existing_budgets[budget_type] = budget
    elseif budget_type == "CVaR"
        if !haskey(model, :u) || !haskey(model, :cvar_term)
            @info("Initializing risk terms for Weighted CVaR objective.")
            add_risk_terms!(model, inputs, settings)
            settings["risk aversion flag"] = true
        end
        @constraint(model, c_budget_cvar, model[:cvar_term] <= budget)
        push!(budget_added, :c_budget_cvar)
        existing_budgets[budget_type] = budget

    elseif budget_type == "System_Expected"
        @constraint(model, c_budget_sys_exp, model[:exp_sys_cost] <= budget)
        push!(budget_added, :c_budget_sys_exp)
        existing_budgets[budget_type] = budget
    elseif budget_type == "System_Weighted_CVaR"
        if risk_aversion < 0 || risk_aversion > 1
            error("For System_CVaR budget constraint, please provide a valid risk aversion weight between 0 and 1.")
        end
        if !haskey(model, :u) || !haskey(model, :cvar_term)
            @info("Initializing risk terms for Weighted CVaR objective.")
            add_risk_terms!(model, inputs, settings)
            settings["risk aversion flag"] = true
        end
        @constraint(model, c_budget_sys_cvar, model[:inv_cost] + (risk_aversion)*model[:expected_alpha] + (1-risk_aversion)*model[:cvar_term]<= budget)
        push!(budget_added, :c_budget_sys_cvar)
        existing_budgets[budget_type] = budget
    else
        error("Unsupported budget type. Please choose from 'Investment', 'Expected', 'CVaR', 'System_Expected', or 'System_CVaR'.")
    end

    @info("Added budget constraint: ", model[budget_added[1]])
    return existing_budgets
end

function manage_cuts(MP::Model, cuts_to_keep::Vector{String})
    # Get all cut constraints in the model
    all_cuts = [con for con in all_constraints(MP, include_variable_in_set_constraints = false) if startswith(string(con), "optimality_cut_") || startswith(string(con), "cvar_tail_cuts_")]
    
    for cut in all_cuts
        if !(name(cut) in cuts_to_keep)
            delete(MP, cut)
        end
    end
    return cuts_to_keep
end

function deactivate_cuts(MP::Model, cuts_to_deactivate::Vector{String})
    for cut_name in cuts_to_deactivate
        set_normalized_rhs(constraint_by_name(MP, cut_name), 1e6) # effectively deactivate the cut by setting a very large RHS
    end
end

function reactivate_cuts(MP::Model, cuts_to_reactivate::Vector{String}, original_rhs::Dict{String, Float64})
    for cut_name in cuts_to_reactivate
        if haskey(MP, Symbol(cut_name)) && haskey(original_rhs, cut_name) && normalized_rhs(MP[Symbol(cut_name)]) >= 1e6
            set_normalized_rhs(MP[Symbol(cut_name)], original_rhs[cut_name]) # reactivate the cut by restoring original RHS
        else
            @warn("Cut ", cut_name, " not found in model or original RHS not recorded. Cannot reactivate.")
        end
    end
end


function fix_capacities!(MP::Model, new_caps::Vector{Float64}, new_line_caps::Vector{Float64})
    for r in 1:length(new_caps)
        fix(MP[:x][r], new_caps[r], force = true)
    end
    for l in 1:length(new_line_caps)
        fix(MP[:x_line][l], new_line_caps[l], force = true)
    end
end

function unfix_capacities!(MP::Model)
    for r in 1:length(MP[:x])
        unfix(MP[:x][r])
    end
    for l in 1:length(MP[:x_line])
        unfix(MP[:x_line][l])
    end
end