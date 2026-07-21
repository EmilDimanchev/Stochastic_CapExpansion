using Test
using Distributed
using Revise, JuMP, Gurobi, HiGHS, Ipopt, DataFrames, CSV, YAML, Random, LinearAlgebra, Combinatorics, Dates, Distributions, Surrogates

include("../src/Stochastic_CapExpansion.jl")
using .Stochastic_CapExpansion

"""
Standalone probe for the capacity reserve margin (CRM) constraint added to the
zonal Benders economic dispatch subproblem (build_subproblem in
src/model/benders/subproblems_economic_dispatch.jl).

Usage:
  julia --project=. tests/cap_reserve_margin_tests.jl

Environment variables:
  SCE_INPUTS_FOLDER   Input folder path (default: inputs/Inputs_30d_1000scen_7tech_2z)
  SCE_SOLVER          Optional solver override, e.g. HiGHS or Gurobi
"""

const DEFAULT_INPUTS = joinpath("inputs", "Inputs_30d_1000scen_7tech_2z")

function build_crm_settings(inputs_folder::String)
    settings = load_settings(inputs_folder)
    settings["CRM flag"] = true
    settings["Parallel flag"] = false
    if haskey(ENV, "SCE_SOLVER")
        settings["Solver"] = ENV["SCE_SOLVER"]
    end
    return settings
end

function solve_with_capacity(inputs, settings, capacity::Vector{Float64}, capacity_line::Vector{Float64})
    SP = build_subproblem(inputs, settings, [1, 1, 1])
    set_parameter_value.(SP[:x], capacity)
    set_parameter_value.(SP[:x_line], capacity_line)
    optimize!(SP)
    @test termination_status(SP) == MOI.OPTIMAL
    return SP
end

function main()
    inputs_folder = get(ENV, "SCE_INPUTS_FOLDER", DEFAULT_INPUTS)
    @info("CRM constraint probe starting")
    @info("Inputs folder: $(inputs_folder)")

    settings = build_crm_settings(inputs_folder)
    inputs = load_input_data(inputs_folder, settings)

    G = inputs["Number of generation resources"]
    O = inputs["Number of storage resources"]
    R = G + O
    L = inputs["Number of lines"]
    gen_names = inputs["Generation resources"]
    derate = inputs["CRM Derating Factors"]

    peak_z1 = maximum(inputs["Demand"][:, 1])
    peak_z2 = maximum(inputs["Demand"][:, 2])
    gas_z1_idx = findfirst(==("Gas-z1"), gen_names)
    gas_z2_idx = findfirst(==("Gas-z2"), gen_names)

    @testset "Capacity Reserve Margin" begin

        @testset "Derating factors load correctly (elementwise, not vector-keyed)" begin
            @test derate isa Dict
            @test derate["Gas-z1"] == 0.93
            @test derate["Wind-z1"] == 0.8
            @test length(derate) == length(gen_names)
        end

        @testset "Ample capacity: margin met without slack" begin
            capacity = zeros(R)
            capacity[gas_z1_idx] = 2.0 * peak_z1
            capacity[gas_z2_idx] = 2.0 * peak_z2
            capacity_line = zeros(L)

            SP = solve_with_capacity(inputs, settings, capacity, capacity_line)

            slack = value.(SP[:vCRMSlack])
            @test all(slack .<= 1e-6)
            crm_dual = dual.(SP[:cCapacityResMargin])
            @test all(crm_dual .>= -1e-6)
        end

        @testset "Tight capacity: margin requires slack, at price cap dual" begin
            # Deliberately starved capacity -- well short of what's needed to
            # clear the 15.6% reserve margin even after voluntary curtailment
            # (segments 2+) is fully exhausted, so the price-cap slack must engage.
            capacity = zeros(R)
            capacity[gas_z1_idx] = 0.5 * peak_z1
            capacity[gas_z2_idx] = 0.5 * peak_z2
            capacity_line = zeros(L)

            SP = solve_with_capacity(inputs, settings, capacity, capacity_line)

            slack = value.(SP[:vCRMSlack])
            @test maximum(slack) > 0.0

            # Where slack is strictly positive, complementary slackness pins the
            # shadow price to exactly the (unscaled) PriceCap for that zone.
            t_weights = inputs["Period weights"][:, 1]
            crm_dual = dual.(SP[:cCapacityResMargin])
            crm_price_cap = inputs["CRM Price Cap"] ./ settings["Scaling factor cost"]
            for t in 1:size(slack, 1), z in 1:size(slack, 2)
                if slack[t, z] > 1e-6
                    @test isapprox(crm_dual[t, z] / t_weights[t], crm_price_cap[z]; atol = 1e-6)
                end
            end
        end

        @testset "CRM flag off: constraint absent, existing behavior unchanged" begin
            settings_off = build_crm_settings(inputs_folder)
            settings_off["CRM flag"] = false
            inputs_off = load_input_data(inputs_folder, settings_off)

            SP = build_subproblem(inputs_off, settings_off, [1, 1, 1])
            @test !haskey(SP, :cCapacityResMargin)
            @test !haskey(SP, :vCRMSlack)
        end
    end
end

main()
