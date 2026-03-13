module Stochastic_CapExpansion

    using JuMP, Gurobi, DataFrames, CSV, Random, LinearAlgebra
    using Revise, YAML

    import JuMP, Gurobi, DataFrames, CSV, YAML, Random, LinearAlgebra

    export build_optimization_model, run_optimization_model, add_budget_constraint!, add_nse_cap!, set_objective!, write_results, write_exploration_results!, load_input_data
    export generate_weights, load_settings

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
    folders = ["inputs", "model", "results", "settings"]
    for folder in folders
        include_all_in_folder(folder)
    end
    cd("..")
end