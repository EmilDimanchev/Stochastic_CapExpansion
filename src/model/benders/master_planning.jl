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

function build_planning_model(inputs, settings; risk_aversion_weight = 0.5)
    # Model
    MP = Model()
    set_silent(MP)

    # Crossover is forced on (both solvers) so the ordinary per-iteration solve always
    # produces a basis. Benders cuts are added as new rows every iteration
    # (add_optimality_cuts!) with no new columns - the textbook dual-simplex warmstart
    # case - but that only pays off if a basis exists from the previous solve for the
    # solver to warm-start from. Crossover=0/run_crossover="off" (the prior setting)
    # produced pure interior-point solutions with no basis at all, defeating any
    # warmstart before it could even be attempted. This does not apply to the
    # regularization solves (level_set_regularization/qtr_regularization), which
    # deliberately use interior-point methods for numerical reasons and temporarily
    # override these settings for their own solve.
    if settings["Solver"] == "HiGHS"
        my_optimizer = optimizer_with_attributes(HiGHS.Optimizer, "solver" => "choose", "run_crossover" => "on", "primal_feasibility_tolerance" => 1e-3, "optimality_tolerance" => 1e-3, "user_bound_scale" => -6)
        set_optimizer(MP, my_optimizer)
    elseif settings["Solver"] == "Gurobi"
        set_optimizer(MP, Gurobi.Optimizer)
        set_optimizer_attribute(MP, "OptimalityTol", 1e-5)
        set_optimizer_attribute(MP, "FeasibilityTol", 1e-3)
        set_optimizer_attribute(MP, "Crossover", 1)
        set_optimizer_attribute(MP, "Method", 1)
    end

    # Sets and flags
    risk_aversion_flag = settings["Risk aversion flag"]
    crm_flag = settings["CRM flag"]

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

    Ω = risk_aversion_weight
    # Per-resource upper bound from the input data (Capacity_UB in Resources.csv/
    # Resources_storage.csv) rather than an arbitrary global 1e6: that constant had no
    # relationship to the data (it's ~5x every resource's own real UB) and inflated the
    # model's RHS/bound coefficient range for no reason. Tightening it to the real,
    # per-resource bound only shrinks that range - it never removes a solution that was
    # reachable before, since 1e6 was never actually binding.
    x_ub = inputs["Capacity upper bounds"]

    # ~~~
    # Model formulation
    # ~~~

    # Capacity
    @variable(MP, x[r in 1:R] >= 1e-5) # Capacity, MW
    @constraint(MP, max_capacity[r in 1:R], x[r] <= x_ub[r])

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
    #risk_aversion_weight = #settings["Risk aversion weight"]

    P = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]

    S = size(P)[1] # number of demand scenarios
    F = size(P_f)[1] # number of gas price scenarios
    K = size(P_k)[1] # number of weather scenarios

    # CVaR parameters
    Ψ = settings["Value-at-Risk percent"] 

    # Auxiliary varliables for CVaR
    @variable(MP, u[s in 1:S, f in 1:F, k in 1:K] >= 0) # loss relative to VaR, $/MW
    @variable(MP, ζ) # VaR variable, $/MW
    
    # CVaR definition constraint
    @constraint(MP, cvar_tail[s in 1:S, f in 1:F, k in 1:K], MP[:u][s,f,k] >= -MP[:ζ])
    
    # Accounting for risk in the objective function
    @expression(MP, cvar_term, ζ + 1/Ψ*sum(P[s]*P_f[f]*P_k[k]*u[s,f,k] for s in 1:S, f in 1:F, k in 1:K))
    #@expression(MP, risk_adjusted_sys_cost, MP[:inv_cost] + Ω*MP[:expected_alpha] + (1-Ω)*MP[:cvar_term])

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

    # Capture the ConstraintRefs JuMP hands back from constraint creation so callers can
    # cache them (name => ConstraintRef) instead of later paying for constraint_by_name,
    # which forces MathOptInterface to rebuild its entire name index (see deactivate_cuts).
    new_cut_refs = Dict{String, ConstraintRef}()
    if cvar_tail
        cvar_cuts = @constraint(MP, [s in 1:S, f in 1:F, k in 1:K], coeffs[s,f,k]*MP[:u][s,f,k] >= SP_obj[s,f,k] + sum(cap_dual[r,s,f,k]*(MP[:x][r]-x_prev[r]) for r in 1:R) + sum(line_dual[l,s,f,k]*(MP[:x_line][l]-x_prev_line[l]) for l in 1:L) - MP[:ζ], base_name = "cvar_tail_cuts_"*case*"_"*string(iteration)*"_")
        for con in cvar_cuts
            new_cut_refs[name(con)] = con
        end
    end

    opt_cuts = @constraint(MP, [s in 1:S, f in 1:F, k in 1:K], coeffs[s,f,k]*MP[:alpha][s,f,k] >= SP_obj[s,f,k] + sum(cap_dual[r,s,f,k]*(MP[:x][r]-x_prev[r]) for r in 1:R) + sum(line_dual[l,s,f,k]*(MP[:x_line][l]-x_prev_line[l]) for l in 1:L), base_name = "optimality_cut_"*case*"_"*string(iteration)*"_")
    for con in opt_cuts
        new_cut_refs[name(con)] = con
    end
    return new_cut_refs
end

#=====================================

Runner

=====================================#

# Multiple dispatch run_planning_model with and without SP outputs for cuts
function run_planning_model(MP, settings, risk_aversion_weight)

    optimize!(MP)
    _log_solver_work!("MP", MP, settings)
    if termination_status(MP) != MOI.OPTIMAL
        @warn("Model did not solve to optimality. Status: ", termination_status(MP))
        unset_silent(MP)
        optimize!(MP)
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
    output = write_outputs(MP, settings, risk_aversion_weight = risk_aversion_weight)

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
    max_gamma = 0.75
    min_gamma = 0.05 # keeps the level-set band from shrinking to a numerically fragile sliver
	if Ia>=0 && Ip>0

		r = Ia/Ip

		if r <=η₁
			#bad iteration: there wasn't enough improvement, we increase γ because we are less confident in the piecewise approximation in the master model
			gamma = gamma/ω;
		elseif η₁ <= r <= η₂
			#okish iteration
			#do nothing
		else
			#good iteration, we are feeling confident in the piecewise approximation in the master model and so we reduce γ
			gamma = gamma*ω;
		end
        gamma = clamp(gamma, min_gamma, max_gamma)

	else
		#do nothing
	end
    return gamma
end


function level_set_regularization(MP, UB, LB, gamma, settings)
    gamma = 0.5
    level_set_rhs = LB + gamma*(UB-LB)
    # Kept as a persistent constraint (RHS updated in place) rather than added and
    # deleted every call: deleting a row forces JuMP/MOI to re-index every constraint
    # after it and discards any basis information for the whole model, whereas
    # set_normalized_rhs is a pure incremental modify. Between regularization calls the
    # bound is relaxed (see below) rather than the row being removed, since the ordinary
    # (non-regularized) master solve that runs in between must not see a stale tight
    # level-set bound from a previous iteration.
    if haskey(MP.obj_dict, :cLevel_set)
        set_normalized_rhs(MP[:cLevel_set], level_set_rhs)
    else
        @constraint(MP, cLevel_set, MP[:eObj] <= level_set_rhs)
    end

    if settings["Solver"] == "HiGHS"
        set_optimizer(MP, () -> Clarabel.Optimizer())
        #set_optimizer(MP, Ipopt.Optimizer)
        #set_optimizer_attribute(MP, "solver", "ipm")
        # Crossover is left off for every other solve in this algorithm (see
        # build_planning_model), so master solves are IPM-only, non-vertex points. That's
        # fine for a plain min-cost solve, but this level-set solve has a degenerate
        # objective (Min 0.0) subject to a tight eObj cap, which is exactly the shape of
        # LP HiGHS's IPM-without-crossover has been observed to spuriously report
        # INFEASIBLE on (confirmed via compute_conflict! finding zero conflicting
        # constraints when this happens - i.e. not a real infeasibility). Turning
        # crossover on for just this solve gets an exact vertex/basic solution instead.
        #set_optimizer_attribute(MP, "run_crossover", "on")
        #set_optimizer_attribute(MP, "presolve", "off")
        set_silent(MP)
    elseif settings["Solver"] == "Gurobi"
        set_optimizer_attribute(MP, "Method", 2)
    end
    #unset_silent(MP)
    @objective(MP, Min, 0.0)
    optimize!(MP)

    if termination_status(MP) != MOI.OPTIMAL && termination_status(MP) != MOI.LOCALLY_SOLVED
        @warn("Regularization problem did not solve to optimality. Status: ", termination_status(MP))
        unset_silent(MP)
        optimize!(MP)
    end
    infeasible = termination_status(MP) == MOI.INFEASIBLE || termination_status(MP) == MOI.INFEASIBLE_OR_UNBOUNDED
    if infeasible
        compute_conflict!(MP)
        list_of_conflicting_constraints = ConstraintRef[];
        for (F, S) in list_of_constraint_types(MP)
            for con in all_constraints(MP, F, S)
                if get_attribute(con, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
                    push!(list_of_conflicting_constraints, con)
                end
            end
        end
        if isempty(list_of_conflicting_constraints)
            @warn("Level-set solve reported infeasible with no conflicting constraints identified - treating as a numerical solver failure rather than a real infeasibility, and falling back to the unstabilized solution for this iteration.")
        else
            display(list_of_conflicting_constraints)
            @warn("Level-set solve is infeasible with identified conflicting constraints (see above); falling back to the unstabilized solution for this iteration.")
        end
    end

    # Write outputs (only meaningful if the level-set solve actually succeeded)
    output = infeasible ? nothing : write_outputs(MP, settings)

    # Relax rather than delete: keeps the row (and MOI's row indexing) stable across
    # calls. Sized relative to UB rather than an arbitrary large constant, to avoid the
    # "huge RHS next to small-scale rows" numerical range issue noted in deactivate_cuts.
    set_normalized_rhs(MP[:cLevel_set], 1e3 * max(abs(UB), 1.0))
    @objective(MP,Min, MP[:eObj])
    if settings["Solver"] == "HiGHS"
        # Restore the same run_crossover="on" base setting build_planning_model uses, so
        # the ordinary solve on the next iteration still gets a basis to warm-start from -
        # otherwise every regularized iteration (i.e. most of them, since this runs
        # whenever gap*100 >= 1) would silently undo that warmstart setup on the way out.
        my_optimizer = optimizer_with_attributes(HiGHS.Optimizer, "solver" => "choose", "run_crossover" => "on")
        set_optimizer(MP, my_optimizer)
    elseif settings["Solver"] == "Gurobi"
        set_optimizer_attribute(MP, "Method", 1)
    end
    set_silent(MP)

    return output

end

# Quadratic trust region regularization
# per https://www.sciencedirect.com/science/article/pii/S0377221719303617
function qtr_regularization(MP, stab_cent_cap, stab_cent_line, gamma_qtr, settings)
    # Sets
    R = length(MP[:x])
    L = length(MP[:x_line])
    @info("Adding trust region constraint with gamma: ", gamma_qtr)
    if settings["Solver"] == "HiGHS"
        set_silent(MP)
        set_optimizer(MP, Ipopt.Optimizer)
    end

    @constraint(MP, cTrust_region, sum((MP[:x][r]-stab_cent_cap[r])^2 for r in 1:R) + sum((MP[:x_line][l]-stab_cent_line[l])^2 for l in 1:L) <= gamma_qtr)
    
    optimize!(MP)

    if termination_status(MP) != MOI.OPTIMAL
        #@warn("Model did not solve to optimality. Status: ", termination_status(MP))
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
    if settings["Solver"] == "HiGHS"
        # Only HiGHS was swapped away (to Ipopt, above - HiGHS has no QCQP support), so
        # only it needs restoring here; the Gurobi branch never left build_planning_model's
        # optimizer instance, and unconditionally rebuilding it (the prior behavior) threw
        # away Gurobi's warm-start state for no reason. Restore the same attributes
        # build_planning_model uses, not a bare HiGHS.Optimizer with no attributes at all.
        my_optimizer = optimizer_with_attributes(HiGHS.Optimizer, "solver" => "choose", "run_crossover" => "on", "primal_feasibility_tolerance" => 1e-3, "optimality_tolerance" => 1e-3, "user_bound_scale" => -6)
        set_optimizer(MP, my_optimizer)
    end
    return output
end


function set_gamma_qtr(gamma_qtr, phi, U_unst, U_prev, x_vector)
    phi_min = 0.1
    kappa = 0.5
    if abs(1- (U_unst/U_prev)) < 0.01
        phi = min(phi_min, kappa*phi)
    end

    gamma_qtr = phi^2 * norm(x_vector, 1)^2

    return gamma_qtr, phi
end




#=====================================

Outputs

=====================================#


function write_outputs(MP, settings; reg=false, risk_aversion_weight=0.5)
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
        output["Risk adjusted system cost"] = (value(MP[:inv_cost]) + risk_aversion_weight*value(MP[:expected_alpha]) + (1 - risk_aversion_weight)*value(MP[:cvar_term]))*scaling_factor_cost
    else
        output["CVaR Loss"] = 0
        output["VaR"] = 0
        output["CVaR term"] = 0
    end
    
    if !reg
        cuts = [con for con in all_constraints(MP, include_variable_in_set_constraints=false) if startswith(name(con), "optimality_cut") || startswith(name(con), "cvar_tail_cuts")]
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





function add_budget_constraint_bendersMP(model::Model, budget::Float64, budget_type::String, existing_budgets::Dict; risk_aversion::Float64 = -1.0, extreme_values::Dict = Dict())
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
    elseif budget_type == "Transformed"
        if !haskey(model, :u) || !haskey(model, :cvar_term)
            error("Cannot run Transformed budget constraint without risk terms. Please initialize risk terms first.")
        end
        transform_cvar = 1/(extreme_values[1.0]["CVaR"] - extreme_values[0.0]["CVaR"])
        transform_sys = 1/(extreme_values[0.0]["System Expected"] - extreme_values[1.0]["System Expected"])
        @constraint(model, c_budget_transformed, transform_sys*(model[:exp_sys_cost]) + transform_cvar*(model[:cvar_term]) <= budget)
        push!(budget_added, :c_budget_transformed)
        tranformed_dict = Dict("budget" => budget, "transform_cvar" => transform_cvar, "transform_sys" => transform_sys)
        existing_budgets[budget_type] = tranformed_dict
    else
        error("Unsupported budget type. Please choose from 'Investment', 'Expected', 'CVaR', 'Transformed','System_Expected', or 'System_CVaR'.")
    end

    @info("Added budget constraint: ", model[budget_added[1]])
    return existing_budgets
end

function update_budget_constraint_bendersMP!(model::Model, new_budget::Float64, budget_type::String, existing_budgets::Dict)
    if budget_type == "Transformed"
        set_normalized_rhs(model[:c_budget_transformed], new_budget)
        existing_budgets["Transformed"]["budget"] = new_budget
    elseif budget_type in ("Investment", "Expected", "CVaR", "System_Expected")
        con_sym = Dict("Investment" => :c_budget_inv, "Expected" => :c_budget_exp, "CVaR" => :c_budget_cvar, "System_Expected" => :c_budget_sys_exp)[budget_type]
        set_normalized_rhs(model[con_sym], new_budget)
        existing_budgets[budget_type] = new_budget
    elseif budget_type == "System_Weighted_CVaR"
        set_normalized_rhs(model[:c_budget_sys_cvar], new_budget)
        existing_budgets[budget_type] = new_budget
    else
        error("Unsupported budget type. Please choose from 'Investment', 'Expected', 'CVaR', 'Transformed','System_Expected', or 'System_CVaR'.")
    end
    @info("Updated budget constraint ($budget_type) RHS to: ", new_budget)
    return existing_budgets
end

function manage_cuts(MP::Model, cuts_to_keep::Vector{String})
    # Get all cut constraints in the model. Filter on name(con), not string(con): the latter
    # pretty-prints the full constraint (every term), which is far more expensive than the
    # O(1) name lookup and runs over every constraint in the model.
    all_cuts = [con for con in all_constraints(MP, include_variable_in_set_constraints = false) if startswith(name(con), "optimality_cut_") || startswith(name(con), "cvar_tail_cuts_")]

    # Batch delete: a scalar delete(MP, cut) per cut forces both HiGHS.jl and Gurobi.jl to
    # redo an O(total constraints in MP) reindex pass - and on Gurobi, a full GRBupdatemodel
    # sync - once per deleted cut, i.e. O(n_deleted * total_constraints) overall. JuMP's
    # vector delete(MP, ::Vector{ConstraintRef}) lets both solvers' MOI wrappers do that
    # reindex/sync exactly once for the whole batch. See deactivate_cuts below for the same
    # pattern with more detail.
    cuts_to_delete = [cut for cut in all_cuts if !(name(cut) in cuts_to_keep)]
    if !isempty(cuts_to_delete)
        delete(MP, cuts_to_delete)
    end
    return cuts_to_keep
end


# Deactivating a cut now means removing it from the model entirely, rather than masking
# it with a huge RHS (-1e6 rows sitting next to ~1e5-scale rows is exactly the kind of
# bad row-bound ratio that was driving level-set regularization to numerical infeasibility).
# Each removed cut's exact constraint (function + set, i.e. its original coefficients and
# RHS) is archived under its name in cut_archive so it can be rebuilt losslessly later via
# reactivate_cuts, instead of being lost.
#
# cut_refs is a name => ConstraintRef cache (populated by add_optimality_cuts! and
# reactivate_cuts) that the caller maintains across iterations. Looking a cut up there is
# O(1); calling constraint_by_name here instead would force MathOptInterface to rebuild its
# entire name index from scratch (it's invalidated on every delete), turning this loop
# quadratic in the number of constraints in MP. Entries are consumed (deleted) as each cut
# is deactivated, so the cache only ever tracks currently-active cuts.
#
# The actual delete is batched (one delete(MP, ::Vector{ConstraintRef}) call after the
# archive loop) rather than called once per cut inside it. Both HiGHS.jl and Gurobi.jl's MOI
# wrappers implement a scalar delete that reindexes every remaining constraint in MP (an
# O(total constraints) walk) on every single call - and Gurobi's scalar path additionally
# forces a GRBupdatemodel sync per call. Deleting N cuts one at a time is therefore
# O(N * total_constraints), which is what made large deactivations slower than solving the
# bigger master problem outright. Their vector delete methods do that reindex/sync exactly
# once for the whole batch, so this is a call-batching change only - it doesn't touch any
# RHS, coefficient, or constraint value, so it carries none of the numerical risk that the
# RHS-loosening approach above ran into.
function deactivate_cuts(MP::Model, cuts_to_deactivate::Vector{String}, cut_archive::Dict{String, Any}, cut_refs::Dict{String, ConstraintRef})
    @info("Deactivating $(length(cuts_to_deactivate)) cuts")
    time_1 = time()
    cons_to_delete = ConstraintRef[]
    names_to_delete = String[]
    for cut_name in cuts_to_deactivate
        con = get(cut_refs, cut_name, nothing)
        if con === nothing
            con = constraint_by_name(MP, cut_name) # fallback for cuts not tracked in cut_refs
        end
        if con === nothing
            continue # already removed (e.g. via manage_cuts)
        end
        cut_archive[cut_name] = constraint_object(con)
        push!(cons_to_delete, con)
        push!(names_to_delete, cut_name)
    end
    if !isempty(cons_to_delete)
        delete(MP, cons_to_delete)
    end
    for cut_name in names_to_delete
        delete!(cut_refs, cut_name)
        #unregister(MP, Symbol(cut_name)) # They aren't registered in the first place bc of anonymous + base name construction! That's a new thing I've learned today - thanks, Claude!
    end
    time_1 = time() - time_1
    @info("Deactivated $(length(cuts_to_deactivate)) cuts in $(round(time_1, digits=2)) seconds")
    return cut_archive
end

# Rebuilds cuts from cut_archive and re-adds them to MP under their original names,
# removing each from the archive as it's restored. Pass `cut_names` to reactivate a
# specific subset; omit it to reactivate everything currently archived. Returns a
# name => ConstraintRef map of everything that was reactivated, so the caller can merge it
# into its cut_refs cache (see deactivate_cuts) instead of relying on constraint_by_name to
# find these constraints again later.
#
# Intended use: at the start of a new problem (e.g. a new risk weight or MGA direction)
# that reuses the same MP, reactivate everything that was pruned while solving previous
# problems so it gets a fresh look. Since reactivated cuts are new to MP again, the normal
# per-iteration cataloguing in benders_algorithm/mga_benders treats them as brand new cuts
# and gives them the usual `cut_deactivation_threshold` iterations of grace before they're
# eligible to be pruned again - so only the ones that still don't bind for the new problem
# get re-removed.
function reactivate_cuts(MP::Model, cut_archive::Dict{String, Any}, cut_names::Union{Vector{String}, Nothing} = nothing)
    @info("Reactivating $(cut_names === nothing ? length(cut_archive) : length(cut_names)) cuts")
    time_1 = time()
    names_to_reactivate = cut_names === nothing ? collect(keys(cut_archive)) : cut_names
    reactivated = Dict{String, ConstraintRef}()
    for cut_name in names_to_reactivate
        if !haskey(cut_archive, cut_name)
            continue
        end
        con = add_constraint(MP, cut_archive[cut_name], cut_name)
        delete!(cut_archive, cut_name)
        reactivated[cut_name] = con
    end
    time_1 = time() - time_1
    @info("Reactivated $(length(reactivated)) cuts in $(round(time_1, digits=2)) seconds")
    return reactivated
end


function fix_capacities!(MP::Model, new_caps::Vector{Float64}, new_line_caps::Vector{Float64})
    for r in 1:length(new_caps)
        fix(MP[:x][r], new_caps[r], force = true)
    end
    for l in 1:length(new_line_caps)
        fix(MP[:x_line][l], new_line_caps[l], force = true)
        #println("Fixing line ", l, " to capacity ", new_line_caps[l])
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

#================

Distributed Sampling Functions

===================#

function _build_mps_for_samples()
    MP = build_planning_model(WORKER_DATA_CACHE[:inputs], WORKER_DATA_CACHE[:settings])
    WORKER_MODEL_CACHE[:MP] = MP
    return nothing
end

function _run_planning_model_for_sample(sample::Vector{Float64})
    MP = WORKER_MODEL_CACHE[:MP]
    sample_caps = sample[1:length(MP[:x])]
    sample_line_caps = sample[length(MP[:x])+1:end]
    fix_capacities!(MP, sample_caps, sample_line_caps)
    output = run_planning_model(MP, WORKER_DATA_CACHE[:settings])
    unfix_capacities!(MP)
    return output
end

function _clear_worker_mp_cache!()
    if WORKER_MODEL_CACHE[:MP] !== nothing
        empty!(WORKER_MODEL_CACHE[:MP])
    end
end

function run_distributed_sampling(samples::Vector{Vector{Float64}})
    # Build model on each worker
    
    @info("Building master problem on each worker...")
    pids = workers()
    n_samples = length(samples)
    samples_per_worker = ceil(Int, n_samples / length(pids))
    if samples_per_worker < 1
        pids = pids[1:n_samples]
    end
    # Build models on workers
    @sync for pid in pids
        @async remotecall_wait(_build_mps_for_samples, pid)
    end
    # Sync and run samples on workers
    outputs = Vector{Dict{String, Any}}(undef, n_samples)
    @sync for (i, sample) in enumerate(samples)
        pid = pids[mod1(i, length(pids))]
        @async begin
            outputs[i] = @fetchfrom pid _run_planning_model_for_sample(sample)
        end
    end
    #@everywhere _clear_worker_mp_cache!()
    @everywhere GC.gc()
    return outputs
end