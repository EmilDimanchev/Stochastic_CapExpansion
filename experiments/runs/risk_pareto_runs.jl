#=================

Risk-pareto, wing-cap and warm-start runs - laptop and Della entry points for run_stochastic_exploration_risk_pareto / run_base_mga

===================#

function risk_pareto_test_della(test_index, budget_type; budget_multiplier = 1.0)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")#joinpath("inputs", "Inputs_30repdays_ext_1000scen_7techs")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)

    configure_parallel_workers!(settings)

    
    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)
    results_folder = joinpath(results_folder, "Pareto5")
    #budget_multiplier = 1.00
    vector_set = nothing
    summary_name = "pareto"
    Eval_SPs = nothing
    mapping=false
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder, summary_folder;
        budget_settings = Dict("Budget type" => budget_type, "Budget multiplier" => budget_multiplier),
        mapping_settings = Dict("Mapping" => false),
        run_settings = Dict("Summary name" => "pareto"))

end

function risk_pareto_test_della_deterministic(test_index, budget_type; budget_multiplier=1.0)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")#joinpath("inputs", "Inputs_30repdays_ext_1000scen_7techs")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)

    configure_parallel_workers!(settings)

    
    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)
    results_folder = joinpath(results_folder, "Pareto5")
    #budget_multiplier = 1.00
    vector_set = nothing
    summary_name = "pareto"
    Eval_SPs = nothing
    mapping=false
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder, summary_folder;
        budget_settings = Dict("Budget type" => budget_type, "Budget multiplier" => budget_multiplier),
        mapping_settings = Dict("Mapping" => false),
        run_settings = Dict("Summary name" => "pareto"))

end

function risk_pareto_test_laptop_5(test_index)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")#joinpath("inputs", "Inputs_30repdays_ext_1000scen_7techs")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)

    configure_parallel_workers!(settings)

    
    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)
    results_folder = joinpath(results_folder, "Pareto5")
    budget_multiplier = 0.00
    vector_set = nothing
    summary_name = "pareto"
    Eval_SPs = nothing
    mapping=true
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder, summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Budget multiplier" => budget_multiplier),
        mapping_settings = Dict("Mapping" => mapping),
        run_settings = Dict("Summary name" => "pareto"))

end

function risk_pareto_tighten_test_laptop(test_index)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)

    configure_parallel_workers!(settings)

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)
    results_folder = joinpath(results_folder, "Pareto5")
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder, summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
        mapping_settings = Dict("Mapping" => true),
        run_settings = Dict("Summary name" => "pareto"))

end

function risk_pareto_tighten_test_della(test_index)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)

    configure_parallel_workers!(settings)

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)
    results_folder = joinpath(results_folder, "Pareto5")
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder, summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
        mapping_settings = Dict("Mapping" => true),
        run_settings = Dict("Summary name" => "pareto"))

end

function risk_pareto_test_no_tighten_della(test_index)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)

    configure_parallel_workers!(settings)

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)
    results_folder = joinpath(results_folder, "Pareto5")
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder, summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => false),
        mapping_settings = Dict("Mapping" => false),
        run_settings = Dict("Summary name" => "pareto"))

end

function wing_caps_test_della(test_index)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z_Della")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)

    configure_parallel_workers!(settings)

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)
    results_folder = joinpath(results_folder, "Pareto5")
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder, summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
        cap_settings = Dict("Cap wings" => true, "Cap levels" => [1, 3]),
        mapping_settings = Dict("Mapping" => true, "Sampling" => false),
        run_settings = Dict("Summary name" => "pareto"))

end



function wing_caps_test_laptop(test_index)

    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)

    configure_parallel_workers!(settings)

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)
    results_folder = joinpath(results_folder, "Pareto5")
    vectors = run_stochastic_exploration_risk_pareto(SPs, inputs, settings, results_folder, summary_folder;
        budget_settings = Dict("Budget type" => "Transformed", "Tighten budget" => true),
        cap_settings = Dict("Cap wings" => true),
        mapping_settings = Dict("Mapping" => false),
        run_settings = Dict("Summary name" => "pareto"))

end


# Minimal diagnostic run for comparing solver warm-start behavior before/after changes to
# Method/Crossover settings (see src/model/benders/master_planning.jl,
# subproblems_economic_dispatch.jl). Skips the full MGA/budget-tightening wrapper - one
# fixed objective, one straight-through benders_algorithm call - so it's fast to re-run
# after each settings change, unlike the full pareto/tighten harness above. Turns on
# "Warmstart debug flag" so each master/subproblem solve logs its simplex/barrier
# iteration count via _log_solver_work! (algorithm.jl) - a more direct signal that warm
# starting is reducing solver work than wall-clock time alone, which is noisy on a laptop.
function warmstart_verification_laptop(test_index)
    inputs_folder = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)
    settings["Warmstart debug flag"] = true

    configure_parallel_workers!(settings)

    MP = build_planning_model(inputs, settings)
    set_objective_bendersMP!(MP, "System_Expected", inputs, settings)
    SPs = build_all_subproblems(inputs, settings)

    elapsed = @elapsed output = benders_algorithm(inputs, settings, MP, SPs, "Warmstart_verification")

    @info("Warmstart verification run complete in $(round(elapsed; digits=2)) seconds")
    @info("Final planning objective: $(output["MP"]["Planning objective"])")
    @info("Final gap: $(output["Gaps"][end])")
    @info("Time MP hist: $(output["Time MP hist"])")
    @info("Time SP hist: $(output["Time SP hist"])")
    @info("Time Reg hist: $(output["Time Reg hist"])")

    return output
end


function simple_comp(test_index)
    inputs_folder = joinpath("inputs","Inputs_30d_1000scen_7tech_2z")
    results_folder = joinpath("outputs", "Test_"*string(test_index))
    summary_folder = joinpath(results_folder, "Summary")
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    if !isdir(summary_folder)
        mkpath(summary_folder)
    end
    settings = load_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)
    
    inputs["Output Demand scenario probabilities"] = inputs["Demand scenario probabilities"] #Establish base weights ahead of time
    inputs["Output Gas price scenario probabilities"] = inputs["Gas price scenario probabilities"]
    inputs["Output Weather scenario probabilities"] = inputs["Weather scenario probabilities"]

    configure_parallel_workers!(settings)

    # Build SPs ------ note that this function set up maintains same SPs across all setups, but each creates its own MP
    SPs = build_all_subproblems(inputs, settings)

    #outputs_mixed, vectors = run_stochastic_exploration(SPs, inputs, settings, joinpath(results_folder, "Both_Flipped"), summary_folder; budget_multiplier=1.10, vector_set=nothing)#vectors
    #outputs_exp, vectors = run_stochastic_exploration_single_type(SPs, inputs, settings, joinpath(results_folder, "Expected"), summary_folder; type = "System_Expected", vector_set = nothing, budget_multiplier=1.10)
    #outputs_cvar, vectors = run_stochastic_exploration_single_type(SPs, inputs, settings, joinpath(results_folder, "CVaR"), summary_folder; type ="System_Weighted_CVaR",vector_set = nothing, budget_multiplier=1.10)
    
    @info("Running base MGA test for mean scenario")
    
    # set all uncertainties to false and risk aversion to false for this test, since we are just running with one scenario selected
    settings["Risk aversion flag"] = false
    settings["Demand uncertainty"] = false
    settings["Gas price uncertainty"] = false
    settings["Weather uncertainty"] = false
    new_inputs = load_input_data(inputs_folder, settings)
    new_inputs["Full Demand scenario probabilities"] = inputs["Demand scenario probabilities"]
    new_inputs["Full Gas price scenario probabilities"] = inputs["Gas price scenario probabilities"]
    new_inputs["Full Weather scenario probabilities"] = inputs["Weather scenario probabilities"]

    i = 0
    results_path = joinpath(results_folder, "Results_Base_MGA", "Scenario_"*string(i))
    _ = run_base_mga(SPs, new_inputs, settings, results_path, summary_folder; budget_multiplier=1.20, vector_set=nothing, scenario=i)

end
