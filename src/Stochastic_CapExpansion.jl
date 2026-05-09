module Stochastic_CapExpansion
    using Revise, JuMP, Gurobi, DataFrames, CSV, YAML, Random, LinearAlgebra, Ipopt, Combinatorics, Dates, Distributed, DistributedArrays, Distributions, Surrogates

   function include_all_in_folder(folder_path::String)
        files = readdir(folder_path, join=true)
        for file in files
            if isdir(file)
                include_all_in_folder(file)
            else
                includet(file)
            end
        end
    end
    cd(@__DIR__)
    println("Loading Stochastic_CapExpansion module from: ", @__DIR__)
    folders = ["inputs", "model", "results", "settings"]
    for folder in folders
        include_all_in_folder(folder)
    end
    cd("..")

    export build_optimization_model, run_optimization_model, add_budget_constraint!, add_nse_cap!, set_objective!, write_results, write_exploration_results!, load_input_data
    export generate_weights, load_settings
    export run_benders_algorithm, benders_algorithm, compute_cvar
    export run_planning_model, build_planning_model, add_optimality_cuts!, add_risk_terms!, add_budget_constraint_bendersMP!, set_objective_bendersMP!
    export run_economic_dispatch, run_all_subproblems, build_all_subproblems, set_capacity_parameters!, run_subproblem
    export write_mapping_results, make_results_mapping_dfs, fix_capacities!, unfix_capacities!
end