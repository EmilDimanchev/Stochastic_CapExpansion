export run_planning_model, build_planning_model, add_optimality_cuts!

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
        output["VaR"] = value.(ζ).*scaling_factor_cost
    else
        output["CVaR Loss"] = zeros(S)
        output["VaR"] = 0
    end

    return output
end


function build_planning_model(inputs, settings)
    # Model
    MP = Model(Gurobi.Optimizer)
    set_silent(MP)
    set_optimizer_attribute(MP, "OptimalityTol", 1e-5)
    set_optimizer_attribute(MP, "FeasibilityTol", 1e-5)
    set_optimizer_attribute(MP, "Crossover", 0)

    # Sets and flags
    risk_aversion_flag = settings["Risk aversion flag"]
    risk_aversion_weight = settings["Risk aversion weight"]

    cost_inv = inputs["Investment costs"]
    P = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]

    scaling_factor_cost = settings["Scaling factor cost"]
    cost_inv = cost_inv./scaling_factor_cost

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

    x_ub = 1e6

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

    # ~~~
    # Model formulation
    # ~~~

    # Capacity
    @variable(MP, x[r in 1:R] >= 1e-5) # Capacity, MW
    @constraint(MP, max_capacity[r in 1:R], x[r] <= x_ub)

    # Cuts
    @variable(MP, alpha[s in 1:S, f in 1:F, k in 1:K] >= 0)

    if risk_aversion_flag
        # Auxiliary varliables for CVaR
        @variable(MP, u[s in 1:S, f in 1:F, k in 1:K] >= 0) # loss relative to VaR, $/MW
        @variable(MP, ζ) # VaR variable, $/MW
    end
        # @constraint(MP, cvar_tail[s in 1:S, f in 1:F, k in 1:K], u[s,f,k] >= sum(t_weights[t]*g[r,t,s,f,k]*cost_var[r,f] for r in 1:G) + sum(t_weights[t]*eNSE_Cost[t,s,f,k] for t in 1:T) + sum(cost_inv[r]*x[r] for r in 1:R) - ζ)
    

    # Cost Expressions
    @expression(MP, inv_cost, sum(x[r]*cost_inv[r] for r in 1:R))
    @expression(MP, expected_alpha, sum(P[s]*P_f[f]*P_k[k]*alpha[s,f,k] for s in 1:S, f in 1:F, k in 1:K))
    if risk_aversion_flag
        @expression(MP, cvar_term, ζ + 1/Ψ*sum(P[s]*P_f[f]*P_k[k]*u[s,f,k] for s in 1:S, f in 1:F, k in 1:K))
    end
    # ~~~
    # Objective function
    # ~~~ 

    if risk_aversion_flag
        @expression(MP, eObj, inv_cost + Ω*expected_alpha + (1-Ω)*cvar_term)
        @expression(MP, eOpObj,Ω*expected_alpha + (1-Ω)*cvar_term) 
        @objective(MP, Min, eObj)
    else
        @expression(MP, eObj, inv_cost + expected_alpha)
        @expression(MP, eOpObj, expected_alpha)
        @objective(MP, Min, eObj)
    end
    return MP

end

function add_optimality_cuts!(MP, SP_obj, SP_dual, x_prev, inputs,settings, iteration)
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

    SP_obj = SP_obj./scaling_factor_cost
    SP_dual = SP_dual./scaling_factor_cost
    if cvar_tail
        @constraint(MP, [s in 1:S, f in 1:F, k in 1:K], MP[:u][s,f,k] >= SP_obj[s,f,k] + sum(SP_dual[r,s,f,k]*(MP[:x][r]-x_prev[r]) for r in 1:R) - MP[:ζ], base_name = "cvar_tail_cuts_"*string(iteration))
    end
    
    @constraint(MP, [s in 1:S, f in 1:F, k in 1:K], MP[:alpha][s,f,k] >= SP_obj[s,f,k] + sum(SP_dual[r,s,f,k]*(MP[:x][r]-x_prev[r]) for r in 1:R), base_name = "optimality_cut_"*string(iteration))
end
# Multiple dispatch run_planning_model with and without SP outputs for cuts
function run_planning_model(MP, settings)

    optimize!(MP)
    println("Planning model objective: ", objective_value(MP))

    # Write outputs
    output = write_outputs(MP, settings)

    return output

end


function regularization(MP, UB, LB, settings)

    @constraint(MP, cLevel_set, MP[:eObj] <= LB + 0.5*(UB-LB))
    
    @objective(MP, Min, sum(0.0 * MP[:alpha]))

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

    return output

end

function write_outputs(MP, settings)
    scaling_factor_cost = settings["Scaling factor cost"]
    output = Dict{String, Any}()
    output["Planning objective"] = value(MP[:eObj])*scaling_factor_cost
    output["Capacity"] = value.(MP[:x])
    output["Alpha"] = value.(MP[:alpha]).*scaling_factor_cost
    output["Contract volume"] = 0
    output["Contract price"] = 0
    output["Inv_cost"] = value(MP[:inv_cost])*scaling_factor_cost
    output["Expected alpha"] = value(MP[:expected_alpha])*scaling_factor_cost
    
    # Record results from the linearization or not
    if settings["Risk aversion flag"]
        output["CVaR Loss"] = value.(MP[:u]).*scaling_factor_cost
        output["VaR"] = value.(MP[:ζ]).*scaling_factor_cost
        output["CVaR term"] = value(MP[:cvar_term])*scaling_factor_cost
    else
        output["CVaR Loss"] = 0
        output["VaR"] = 0
        output["CVaR term"] = 0
    end

    return output

end