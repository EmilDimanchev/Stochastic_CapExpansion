using JuMP
using CSV
using DataFrames
using Gurobi

function run_planning_model(inputs, settings, SP_obj_list, SP_dual_list, x_prev)

    # Model
    MP = Model(Gurobi.Optimizer)
    # set_optimizer_attribute(MP, "OutputFlag", 0)

    # ~~~~
    # Load inputs
    # ~~~

    SP_obj = reduce((a, b) -> cat(a, b; dims=4), SP_obj_list)
    SP_dual = reduce((a, b) -> cat(a, b; dims=5), SP_dual_list)

    risk_aversion_flag = settings["Risk aversion flag"]
    risk_aversion_weight = settings["Risk aversion weight"]

    cost_inv = inputs["Investment costs"]
    P = inputs["Demand scenario probabilities"]
    P_f = inputs["Gas price scenario probabilities"]
    P_k = inputs["Weather scenario probabilities"]

    scaling_factor_cost = settings["Scaling factor cost"]
    cost_inv = cost_inv./scaling_factor_cost
    SP_obj = SP_obj./scaling_factor_cost
    SP_dual = SP_dual./scaling_factor_cost
    x_ub = 1e6

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

    # ~~~
    # Model formulation
    # ~~~

    # Capacity
    @variable(MP, x[r in 1:R] >= 0) # Capacity, MW
    @constraint(MP, max_capacity[r in 1:R], x[r] <= x_ub)

    if risk_aversion_flag
        # Auxiliary varliables for CVaR
        @variable(MP, u[s in 1:S, f in 1:F, k in 1:K] >= 0) # loss relative to VaR, $/MW
        @variable(MP, ζ) # VaR variable, $/MW
        # @constraint(MP, cvar_tail[s in 1:S, f in 1:F, k in 1:K], u[s,f,k] >= sum(t_weights[t]*g[r,t,s,f,k]*cost_var[r,f] for r in 1:G) + sum(t_weights[t]*eNSE_Cost[t,s,f,k] for t in 1:T) + sum(cost_inv[r]*x[r] for r in 1:R) - ζ)
        @constraint(MP, cvar_tail[s in 1:S, f in 1:F, k in 1:K, j in 1:J], u[s,f,k] >= SP_obj[s,f,k,j] + sum(SP_dual[r,s,f,k,j]*(x[r]-x_prev[j][r]) for r in 1:R) - ζ)
    end 

    # Cuts
    @variable(MP, alpha[s in 1:S, f in 1:F, k in 1:K] >= 0)

    @constraint(MP, optimality_cuts[s in 1:S, f in 1:F, k in 1:K, j in 1:J], alpha[s,f,k] >= SP_obj[s,f,k,j] + sum(SP_dual[r,s,f,k,j]*(x[r]-x_prev[j][r]) for r in 1:R))

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