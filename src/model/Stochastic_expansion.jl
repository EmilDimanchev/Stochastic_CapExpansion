using JuMP
using CSV
using DataFrames
using Gurobi
using Random
using LinearAlgebra

function set_objective!(model, objective_type::String; obj_weight::Float64 = -1.0, set_coeffs::Vector = [])
    if objective_type == "System_Expected"
        @objective(model, Min, model[:investment_cost] + model[:expected_cost])
        @info("Objective set to minimize expected system cost: ", objective_function(model))
    elseif objective_type == "System_Weighted_CVaR"
        if obj_weight < -0.0001 || obj_weight > 1.001
            error("For Weighted CVaR objective, please provide a valid weight between 0 and 1.")
        end
        @objective(model, Min, model[:investment_cost] + obj_weight*model[:expected_cost] + (1-obj_weight)*model[:cvar])
        @info("Objective set to minimize weighted CVaR. Weight is: ", obj_weight, " Objective function: ", objective_function(model))
    elseif objective_type == "Capacity"
        if set_coeffs == []
            error("For Capacity objective, please provide a vector of coefficients for the capacity expression.")
        end
        @objective(model, Min, set_coeffs'*model[:x])
        @info("Objective set to minimize capacity expression: ", objective_function(model))
    else
        error("Unsupported objective type. Please choose from 'Expectation', 'Weighted CVaR', or 'Capacity'.")
    end
end

function generate_weights(iterations::Int64, num_resources::Int, option::String, settings::Dict)
    if settings["Seed"] !== nothing
        Random.seed!(settings["Seed"])
    end
    weights = Vector{Vector{Float64}}(undef, 0)
    vector = zeros(num_resources)
    for i in 1:ceil(Int, iterations/2)
        if option == "random"
            vector = randn(num_resources)*1000
        elseif option == "min/max"
            vector = rand(-1:1, num_resources)*1000
        end
        push!(weights, vector)
        push!(weights, -1*vector) 
    end

    weights = sort_weights(weights)

    return weights
end

function sort_weights(weights::Vector{Vector{Float64}})
    norms = [norm(weights[i]) for i in 1:length(weights)]
    dot_product = collect(dot(weights[i], weights[1]) for i in 1:length(weights))./norms
    dot_product = clamp.(dot_product, -1.0, 1.0)
    angles = [acos(dot_product[i]) for i in 1:length(weights)]
    sorted_indices = sortperm(angles)
    weights = weights[sorted_indices]

    return weights
end

# Add budget constraint

function add_budget_constraint!(model::Model, budget::Float64, budget_type::String; risk_aversion::Float64 = -1.0)
    budget_added = Symbol[]
    if budget_type == "Investment"
        @constraint(model, c_budget_inv, model[:investment_cost] <= budget)
        push!(budget_added, :c_budget_inv)
    elseif budget_type == "Expected"
        @constraint(model, c_budget_exp, model[:expected_cost] <= budget)
        push!(budget_added, :c_budget_exp)
    elseif budget_type == "CVaR"
        @constraint(model, c_budget_cvar, model[:cvar] <= budget)
        push!(budget_added, :c_budget_cvar)
    elseif budget_type == "System_Expected"
        @constraint(model, c_budget_sys_exp, model[:investment_cost] + model[:expected_cost] <= budget)
        push!(budget_added, :c_budget_sys_exp)
    elseif budget_type == "System_Weighted_CVaR"
        if risk_aversion < 0 || risk_aversion > 1
            error("For System_CVaR budget constraint, please provide a valid risk aversion weight between 0 and 1.")
        end
        @constraint(model, c_budget_sys_cvar, model[:investment_cost] + (risk_aversion)*model[:expected_cost] + (1-risk_aversion)*model[:cvar]<= budget)
        push!(budget_added, :c_budget_sys_cvar)
    else
        error("Unsupported budget type. Please choose from 'Investment', 'Expected', 'CVaR', 'System_Expected', or 'System_CVaR'.")
    end

    @info("Added budget constraint: ", model[budget_added[1]])
    return model
end

function add_budget_constraint!(model::Model, budget::Vector, budget_type::Vector; risk_aversion::Float64 = -1.0)
    budget_added = Symbol[]
    for i in 1:length(budget)
        if budget_type[i] == "Investment"
            @constraint(model, c_budget_inv, model[:investment_cost] <= budget[i])
            push!(budget_added, :c_budget_inv)
        elseif budget_type[i] == "Expected"
            @constraint(model, c_budget_exp, model[:expected_cost] <= budget[i])
            push!(budget_added, :c_budget_exp)
        elseif budget_type[i] == "CVaR"
            @constraint(model, c_budget_cvar, model[:cvar] <= budget[i])
            push!(budget_added, :c_budget_cvar)
        elseif budget_type[i] == "System_Expected"
            @constraint(model, c_budget_sys_exp, model[:investment_cost] + model[:expected_cost] <= budget[i])
            push!(budget_added, :c_budget_sys_exp)
        elseif budget_type[i] == "System_Weighted_CVaR"
            if risk_aversion < 0 || risk_aversion > 1
                error("For System_CVaR budget constraint, please provide a valid risk aversion weight between 0 and 1.")
            end
            @constraint(model, c_budget_sys_cvar, model[:investment_cost] + (1-risk_aversion)*model[:cvar] + (risk_aversion)*model[:expected_cost] <= budget[i])
            push!(budget_added, :c_budget_sys_cvar)
        else
            error("Unsupported budget type. Please choose from 'Investment', 'Expected', 'CVaR', 'System_Expected', or 'System_CVaR'.")
        end

        @info("Added budget constraint: ", model[budget_added[i]])
    end
    return model
end

function add_nse_cap!(model::Model, cap::Float64)
    @constraint(model, cNSE_Cap, sum(model[:eNSE][t,s,f] for t in 1:size(model[:eNSE],1), s in 1:size(model[:eNSE],2), f in 1:size(model[:eNSE],3)) <= cap)
end