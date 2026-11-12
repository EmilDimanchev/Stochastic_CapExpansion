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

    # Solver attributes live in master_highs_optimizer / master_gurobi_attributes, so the
    # regularization solves (which use interior-point methods and override some of them for
    # their own solve) restore exactly what is set here.
    if settings["Solver"] == "HiGHS"
        set_optimizer(MP, master_highs_optimizer())
    elseif settings["Solver"] == "Gurobi"
        set_optimizer(MP, Gurobi.Optimizer)
        set_gurobi_attributes!(MP, master_gurobi_attributes(settings))
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
    # A resource whose lower bound meets its upper bound (a fleet that can't retire) is fixed rather
    # than boxed: a box with no interior is a problem for the barrier regularization solves.
    x_lb = inputs["Capacity lower bounds"]
    fixed_resources = x_lb .>= x_ub
    @variable(MP, x[r in 1:R] >= max(1e-5, x_lb[r])) # Capacity, MW
    fix.(x[fixed_resources], x_ub[fixed_resources]; force = true)
    @constraint(MP, max_capacity[r in findall(.!fixed_resources)], x[r] <= x_ub[r])

    @variable(MP, x_line[l in 1:L] >= 0) # Transmission line expansion, MW
    @constraint(MP, max_line_expansion[l in 1:L], x_line[l] <= inputs["Max expansion"][l])

    # Operating cost estimates: one per subproblem, i.e. per scenario and time block (see
    # subproblem_blocks). Each block gets its own optimality cuts (add_optimality_cuts!), and a
    # scenario's estimate is the sum over its blocks. docs/benders_block_cuts.pdf has the math.
    B = length(subproblem_blocks(inputs, settings))
    @variable(MP, alpha[s in 1:S, f in 1:F, k in 1:K, b in 1:B] >= 0)
    @expression(MP, scenario_alpha[s in 1:S, f in 1:F, k in 1:K], sum(alpha[s,f,k,b] for b in 1:B))


    # Cost Expressions
    @expression(MP, inv_cost, sum(x[r]*cost_inv[r] for r in 1:R) + sum(x_line[l]*capex_line[l] for l in 1:L))
    @expression(MP, expected_alpha, sum(P[s]*P_f[f]*P_k[k]*scenario_alpha[s,f,k] for s in 1:S, f in 1:F, k in 1:K))
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
        MP.ext[:objective_risk_weight] = Float64(Ω)
    else
        @expression(MP, eObj, inv_cost + expected_alpha)
        @objective(MP, Min, eObj)
        MP.ext[:objective_risk_weight] = 1.0
    end
    return MP

end

# Weight w of the master's cost objective, inv + w·E[op cost] + (1 - w)·CVaR[op cost] (1.0 for expected
# cost alone). build_planning_model and set_objective_bendersMP! record it whenever they set eObj, and
# benders_algorithm builds LB and UB from it, so both bounds always describe the objective the master
# actually minimizes. Not a settings flag: "Risk aversion flag" can be false while the objective has a
# CVaR term (set_objective_bendersMP! adds the risk terms itself), and bounds built from the flag then
# compare inv + E[cost] with a master minimizing the weighted objective, which goes negative.
# nothing: the objective isn't a cost (MGA capacity directions).
objective_risk_weight(MP::Model) = get(MP.ext, :objective_risk_weight, nothing)

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
    
    # CVaR definition constraint: u is each scenario's cost estimate above VaR. Tied to the block cost
    # estimates (scenario_alpha), so CVaR sees every block cut and no CVaR cut rows are needed. With
    # alpha >= 0 this also covers the old u >= -ζ row. See docs/benders_block_cuts.pdf.
    @constraint(MP, cvar_tail[s in 1:S, f in 1:F, k in 1:K], MP[:u][s,f,k] >= MP[:scenario_alpha][s,f,k] - MP[:ζ])
    
    # Accounting for risk in the objective function
    @expression(MP, cvar_term, ζ + 1/Ψ*sum(P[s]*P_f[f]*P_k[k]*u[s,f,k] for s in 1:S, f in 1:F, k in 1:K))
    #@expression(MP, risk_adjusted_sys_cost, MP[:inv_cost] + Ω*MP[:expected_alpha] + (1-Ω)*MP[:cvar_term])

end

# One optimality cut per subproblem (scenario × time block) at the capacity point (x_prev, x_prev_line):
#     coeff·alpha[s,f,k,b] >= θ_b + β_bᵀ(x − x_prev) + β_line_bᵀ(x_line − x_prev_line)
# SP_obj and coeffs are S × F × K × B, cap_dual R × S × F × K × B, line_dual L × S × F × K × B, all
# unscaled (the "Block ..." terms of run_all_subproblems). A block that didn't solve has coeff 0, so
# its row carries no information. CVaR needs no rows of its own: cvar_tail ties u to the alphas.
function add_optimality_cuts!(MP, SP_obj, cap_dual, line_dual, x_prev, x_prev_line, coeffs, inputs, settings, iteration, case)
    scaling_factor_cost = settings["Scaling factor cost"]
    S, F, K, B = size(MP[:alpha])
    size(SP_obj) == (S, F, K, B) || error("add_optimality_cuts!: SP_obj is $(size(SP_obj)), the master's alpha is $((S, F, K, B)).")
    R = length(MP[:x])
    L = length(MP[:x_line])

    SP_obj = SP_obj./scaling_factor_cost
    cap_dual = cap_dual./scaling_factor_cost
    line_dual = line_dual./scaling_factor_cost

    # Capture the ConstraintRefs JuMP hands back from constraint creation so callers can
    # cache them (name => ConstraintRef) instead of later paying for constraint_by_name,
    # which forces MathOptInterface to rebuild its entire name index (see CutPool).
    new_cut_refs = Dict{String, ConstraintRef}()
    opt_cuts = @constraint(MP, [s in 1:S, f in 1:F, k in 1:K, b in 1:B], coeffs[s,f,k,b]*MP[:alpha][s,f,k,b] >= SP_obj[s,f,k,b] + sum(cap_dual[r,s,f,k,b]*(MP[:x][r]-x_prev[r]) for r in 1:R) + sum(line_dual[l,s,f,k,b]*(MP[:x_line][l]-x_prev_line[l]) for l in 1:L), base_name = "optimality_cut_"*case*"_"*string(iteration)*"_")
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
    restore_master_basis!(MP) # no-op unless a level-set pass saved one
    optimize!(MP)
    _log_solver_work!("MP", MP, settings)
    if termination_status(MP) != MOI.OPTIMAL
        # One verbose retry. Its result is used if it solves: dropping it used to leave an empty
        # output and a KeyError downstream ("Expected alpha", 16-zone run, Nov 2026).
        @warn("Model did not solve to optimality. Status: ", termination_status(MP))
        unset_silent(MP)
        optimize!(MP)
        set_silent(MP)
    end
    if termination_status(MP) == MOI.OPTIMAL
        output = write_outputs(MP, settings, risk_aversion_weight = risk_aversion_weight)
    end
    if termination_status(MP) == MOI.INFEASIBLE_OR_UNBOUNDED || termination_status(MP) == MOI.INFEASIBLE
        @warn("Model did not solve to optimality. Status: ", termination_status(MP), "Conflicting constraints: ")
        compute_conflict!(MP)
        list_of_conflicting_constraints = ConstraintRef[];
        for (F, S) in list_of_constraint_types(MP)
            for con in all_constraints(MP, F, S)
                if get_attribute(con, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
                    push!(list_of_conflicting_constraints, con)
                end
            end
        end
        show_conflicts(list_of_conflicting_constraints)
        error("Model infeasible. See conflicting constraints above.")
    end
    termination_status(MP) == MOI.OPTIMAL || error("Master problem did not solve to optimality, even on the retry. Status: $(termination_status(MP))")

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


# The HiGHS options every ordinary master solve uses. The level-set and QTR solves swap the
# optimizer out and back, and set_optimizer drops every attribute, so they restore from here.
master_highs_optimizer() = optimizer_with_attributes(HiGHS.Optimizer, "solver" => "choose", "run_crossover" => "off",
    "primal_feasibility_tolerance" => 1e-3, "optimality_tolerance" => 1e-3, "user_bound_scale" => -6)

# The Gurobi attributes of every ordinary master solve. level_set_regularization switches to
# level_set_gurobi_attributes for its own solve and back to these afterwards, so both the build and
# the restore read this one list.
# Concurrent (Method 3) so the master uses every core: dual and primal simplex warm-start from the
# previous basis (LPWarmStart 1 hands the basis to simplex only) while barrier runs alongside, and the
# first to finish wins. Dual simplex alone (the old Method 1) is single-threaded, and the master was
# 93-96% of wall time on the 16-zone case with the other cores idle. Crossover stays on (-1) so a
# barrier win still leaves a basis for the next iteration's warm start.
# "Threads" in "Master problem settings" caps the threads: Gurobi's default (0) counts every logical core on the node, not
# the cores Slurm allocated to the job.
master_gurobi_attributes(settings) = [
    "OptimalityTol" => 1e-5, "FeasibilityTol" => 1e-3, "MultiObjMethod" => 0,
    "Method" => 3, "LPWarmStart" => 1, "Crossover" => -1, "Threads" => master_settings(settings)["Threads"]]

# Level-set solve: barrier for an interior point, so crossover is off. LPWarmStart 0, since with a
# basis from the previous solve LPWarmStart 2 makes Gurobi skip barrier even with Method 2 (see
# level_set_regularization).
level_set_gurobi_attributes(settings) = [
    "Method" => 2, "LPWarmStart" => 0, "Crossover" => 0, "Threads" => master_settings(settings)["Threads"]]

set_gurobi_attributes!(model, attributes) = foreach(((k, v),) -> set_optimizer_attribute(model, k, v), attributes)

# Adds the level-set row eObj <= rhs. level_set_regularization deletes it again after its solve, so
# each call builds it from the current eObj: one MP is reused across the risk-weight anchors, and
# set_objective_bendersMP! replaces eObj with a new expression each time.
# Deleted, not relaxed: an Inf RHS reaches Gurobi as 1e100, which Gurobi treats as a real (huge) RHS.
# On the 16-zone case that row sent the concurrent master's barrier to an objective of -5e99 and a
# DUAL_INFEASIBLE status. A large finite RHS (1e3·|UB|) widened the RHS range the same way.
function add_level_set_row!(MP, rhs)
    haskey(MP.obj_dict, :cLevel_set) && error("add_level_set_row!: the master already has a level-set row.")
    return @constraint(MP, cLevel_set, MP[:eObj] <= rhs)
end

function delete_level_set_row!(MP)
    delete(MP, MP[:cLevel_set])
    unregister(MP, :cLevel_set)
end

#=====================================
Master basis across a level-set solve (Gurobi)

The level-set solve is barrier with crossover off, so it leaves no basis, and the next master
solve would start simplex cold. save_master_basis! records the basis of the last master solve
before the level-set solve, and run_planning_model hands it back to Gurobi (restore_master_basis!)
just before the next master solve. Rows added in between (new cuts, separation restores) get a
basic slack, and rows deleted in between (pruned cuts, the level-set row) are dropped: the same
start the master gets on an iteration without regularization. If Gurobi finds the basis
inconsistent it ignores it and starts cold, which is no worse than not saving it.

HiGHS doesn't need this: its level-set solve swaps the optimizer with set_optimizer, which
rebuilds the solver model, so there is no basis left to keep.
=====================================#

const _VBASIS = "VBasis"
const _CBASIS = "CBasis"
const _BASIC = 0 # Gurobi's VBasis/CBasis code for basic

# Linear rows of the master; these are the ones Gurobi keeps a CBasis for.
master_rows(MP) = [con for (F, S) in list_of_constraint_types(MP) if F <: GenericAffExpr
                   for con in all_constraints(MP, F, S)]

# Must run right after a master solve, before any row is added or deleted: Gurobi discards the
# basis once the model changes.
function save_master_basis!(MP, settings)
    settings["Solver"] == "Gurobi" || return nothing
    vbasis = Dict(x => MOI.get(MP, Gurobi.VariableAttribute(_VBASIS), x) for x in all_variables(MP))
    cbasis = Dict(con => MOI.get(MP, Gurobi.ConstraintAttribute(_CBASIS), con) for con in master_rows(MP))
    MP.ext[:saved_basis] = (vbasis = vbasis, cbasis = cbasis)
    return nothing
end

# The CBasis to set on each row now in the master: its saved status, or basic for a row added
# since the save. Saved rows that have since been deleted aren't in `rows`, so they drop out.
row_basis_statuses(saved_cbasis, rows) = [con => get(saved_cbasis, con, _BASIC) for con in rows]

# Sets a basis saved by save_master_basis! as the start of the next solve, then forgets it.
# A no-op when nothing was saved.
function restore_master_basis!(MP)
    saved = pop!(MP.ext, :saved_basis, nothing)
    saved === nothing && return nothing
    for (x, status) in saved.vbasis
        is_valid(MP, x) && MOI.set(MP, Gurobi.VariableAttribute(_VBASIS), x, status)
    end
    for (con, status) in row_basis_statuses(saved.cbasis, master_rows(MP))
        MOI.set(MP, Gurobi.ConstraintAttribute(_CBASIS), con, status)
    end
    return nothing
end

function level_set_regularization(MP, UB, LB, gamma, settings)
    gamma = 0.5
    add_level_set_row!(MP, LB + gamma*(UB-LB))

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
        # Barrier, without the warm start. With LPWarmStart = 2 and a basis left by the previous
        # master solve, Gurobi skips barrier entirely, even with Method = 2: it crushes that basis
        # and goes straight to crossover and primal simplex. On the 16-zone case that start had
        # primal infeasibility ~7e9, simplex reached Kappa ~1e16, and the solve was declared
        # infeasible. Crossover is off: this solve is meant to return an interior point. All of it is
        # restored after the solve.
        set_gurobi_attributes!(MP, level_set_gurobi_attributes(settings))
    end
    #unset_silent(MP)
    @objective(MP, Min, 0.0)
    log_phase("level-set solve")
    optimize!(MP)

    if termination_status(MP) != MOI.OPTIMAL && termination_status(MP) != MOI.LOCALLY_SOLVED
        @warn("Regularization problem did not solve to optimality. Status: ", termination_status(MP))
        unset_silent(MP)
        optimize!(MP)
    end
    infeasible = termination_status(MP) == MOI.INFEASIBLE 
    numerical_trouble = termination_status(MP) == MOI.INFEASIBLE_OR_UNBOUNDED || termination_status(MP) == MOI.NUMERICAL_ERROR
    if numerical_trouble
        @warn("Numerical trouble encountered, re-running with numerically stable solver")
        if settings["Solver"] == "Gurobi"
            set_optimizer_attribute(MP, "BarHomogeneous", 1)
        end
        optimize!(MP)
    end
    if infeasible
        list_of_conflicting_constraints = ConstraintRef[];
        # The level-set solve runs on Clarabel when Solver == "HiGHS", which doesn't implement
        # compute_conflict!. Treat that as "no conflict identified" (the fallback branch below)
        # instead of crashing the whole run on a diagnostic.
        try
            compute_conflict!(MP)
            for (F, S) in list_of_constraint_types(MP)
                for con in all_constraints(MP, F, S)
                    if get_attribute(con, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
                        push!(list_of_conflicting_constraints, con)
                    end
                end
            end
        catch e
            e isa ArgumentError || rethrow()
            @warn("Conflict analysis unavailable for the current optimizer: $(sprint(showerror, e))")
        end
        if isempty(list_of_conflicting_constraints)
            @warn("Level-set solve reported infeasible with no conflicting constraints identified - treating as a numerical solver failure rather than a real infeasibility, and falling back to the unstabilized solution for this iteration.")
        else
            show_conflicts(list_of_conflicting_constraints)
            @warn("Level-set solve is infeasible with identified conflicting constraints (see above); falling back to the unstabilized solution for this iteration.")
        end
    end

    # Write outputs (only meaningful if the level-set solve actually succeeded)
    output = (infeasible || numerical_trouble) ? nothing : write_outputs(MP, settings)

    # Deleted rather than relaxed to +Inf (see add_level_set_row!)
    delete_level_set_row!(MP)
    @objective(MP,Min, MP[:eObj])
    if settings["Solver"] == "HiGHS"
        # Back to the ordinary master settings. A bare HiGHS.Optimizer here (as before) silently
        # dropped the master's tolerances and user_bound_scale after the first regularized iteration.
        set_optimizer(MP, master_highs_optimizer())
    elseif settings["Solver"] == "Gurobi"
        set_gurobi_attributes!(MP, master_gurobi_attributes(settings))
        if numerical_trouble
            set_optimizer_attribute(MP, "BarHomogeneous", 0)
        end
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
    log_detail("Adding trust region constraint with gamma: $gamma_qtr")
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
        show_conflicts(list_of_conflicting_constraints)
        error("Model infeasible. See conflicting constraints above.")
    end
    
    # Write outputs
    output = write_outputs(MP, settings)

    delete(MP,MP[:cTrust_region])
    unregister(MP,:cTrust_region)
    if settings["Solver"] == "HiGHS"
        # Only HiGHS was swapped away (to Ipopt, above - HiGHS has no QCQP support), so
        # only it needs restoring here; the Gurobi branch never left build_planning_model's
        # optimizer instance, and unconditionally rebuilding it (the prior behavior) threw
        # away Gurobi's warm-start state for no reason. Restore the same attributes
        # build_planning_model uses, not a bare HiGHS.Optimizer with no attributes at all.
        set_optimizer(MP, master_highs_optimizer())
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


function write_outputs(MP, settings; risk_aversion_weight=0.5)
    scaling_factor_cost = settings["Scaling factor cost"]
    output = Dict{String, Any}()
    output["Planning objective"] = value(MP[:eObj])*scaling_factor_cost
    output["Capacity"] = value.(MP[:x])
    output["Line expansion"] = value.(MP[:x_line])
    output["Alpha"] = value.(MP[:scenario_alpha]).*scaling_factor_cost   # per scenario (sum over blocks)
    output["Block alpha"] = value.(MP[:alpha]).*scaling_factor_cost      # S × F × K × B
    output["Contract volume"] = 0
    output["Contract price"] = 0
    output["Inv_cost"] = value(MP[:inv_cost])*scaling_factor_cost
    output["Expected alpha"] = value(MP[:expected_alpha])*scaling_factor_cost
    output["Expected System Cost"] = value.(MP[:exp_sys_cost])*scaling_factor_cost
    output["Inv cost by zone"] = value.(MP[:inv_cost_by_zone]).*scaling_factor_cost
    
    # CVaR results whenever the master has the CVaR terms, which set_objective_bendersMP! and the budget
    # constraints can add even when "Risk aversion flag" is false
    if haskey(MP, :cvar_term)
        output["CVaR Loss"] = value.(MP[:u]).*scaling_factor_cost
        output["VaR"] = value.(MP[:ζ]).*scaling_factor_cost
        output["CVaR term"] = value(MP[:cvar_term])*scaling_factor_cost
        output["Risk adjusted system cost"] = (value(MP[:inv_cost]) + risk_aversion_weight*value(MP[:expected_alpha]) + (1 - risk_aversion_weight)*value(MP[:cvar_term]))*scaling_factor_cost
    else
        output["CVaR Loss"] = 0
        output["VaR"] = 0
        output["CVaR term"] = 0
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
        model.ext[:objective_risk_weight] = 1.0
        log_detail("Objective set to minimize expected system cost: $(objective_function(model))")
    elseif objective_type == "System_Weighted_CVaR"
        if obj_weight < -0.0001 || obj_weight > 1.001
            error("For Weighted CVaR objective, please provide a valid weight between 0 and 1.")
        end
        # Check if model is initialized with risk terms?
        if !haskey(model, :u) || !haskey(model, :cvar_term)
            log_detail("Initializing risk terms for Weighted CVaR objective.")
            add_risk_terms!(model, inputs, settings)
        end
        unregister(model, :eObj)
        @expression(model, eObj, model[:inv_cost] + (obj_weight)*model[:expected_alpha] + (1-obj_weight)*model[:cvar_term])
        @objective(model, Min, model[:eObj])
        model.ext[:objective_risk_weight] = Float64(obj_weight)
        log_detail("Objective set to minimize weighted CVaR. Weight is: $obj_weight Objective function: $(objective_function(model))")
    elseif objective_type == "Capacity"
        if set_coeffs == []
            error("For Capacity objective, please provide a vector of coefficients for the capacity expression.")
        end
        x_full = [collect(model[:x][r] for r in 1:length(model[:x])); collect(model[:x_line][l] for l in 1:length(model[:x_line]))]
        @objective(model, Min, set_coeffs' * x_full)
        delete!(model.ext, :objective_risk_weight)
        log_detail("Objective set to minimize capacity expression: $(objective_function(model))")
    else
        error("Unsupported objective type. Please choose from 'Expectation', 'Weighted CVaR', or 'Capacity'.")
    end
end





function add_budget_constraint_bendersMP(model::Model, budget::Float64, budget_type::String, existing_budgets::Dict; risk_aversion::Float64 = -1.0, extreme_values::Dict = Dict(), constraint_dict::Dict = Dict())
    budget_added = Symbol[]
    if budget_type == "Investment"
        new_budget = @constraint(model, c_budget_inv, model[:inv_cost] <= budget)
        push!(budget_added, :c_budget_inv)
        existing_budgets[budget_type] = budget
        constraint_dict[budget_type] = new_budget
    elseif budget_type == "Expected"
        new_budget = @constraint(model, c_budget_exp, model[:expected_alpha] <= budget)
        push!(budget_added, :c_budget_exp)
        existing_budgets[budget_type] = budget
        constraint_dict[budget_type] = new_budget
    elseif budget_type == "CVaR"
        if !haskey(model, :u) || !haskey(model, :cvar_term)
            log_detail("Initializing risk terms for Weighted CVaR objective.")
            add_risk_terms!(model, inputs, settings)
        end
        new_budget = @constraint(model, c_budget_cvar, model[:cvar_term] <= budget)
        push!(budget_added, :c_budget_cvar)
        existing_budgets[budget_type] = budget
        constraint_dict[budget_type] = new_budget
    elseif budget_type == "System_Expected"
        new_budget = @constraint(model, c_budget_sys_exp, model[:exp_sys_cost] <= budget)
        push!(budget_added, :c_budget_sys_exp)
        existing_budgets[budget_type] = budget
        constraint_dict[budget_type] = new_budget
    elseif budget_type == "System_Weighted_CVaR"
        if risk_aversion < 0 || risk_aversion > 1
            error("For System_CVaR budget constraint, please provide a valid risk aversion weight between 0 and 1.")
        end
        if !haskey(model, :u) || !haskey(model, :cvar_term)
            log_detail("Initializing risk terms for Weighted CVaR objective.")
            add_risk_terms!(model, inputs, settings)
        end
        new_budget = @constraint(model, c_budget_sys_cvar, model[:inv_cost] + (risk_aversion)*model[:expected_alpha] + (1-risk_aversion)*model[:cvar_term]<= budget)
        push!(budget_added, :c_budget_sys_cvar)
        existing_budgets[budget_type] = budget
        constraint_dict[budget_type] = new_budget
    elseif budget_type == "Transformed"
        if !haskey(model, :u) || !haskey(model, :cvar_term)
            error("Cannot run Transformed budget constraint without risk terms. Please initialize risk terms first.")
        end
        transform_cvar = 1/(extreme_values[1.0]["CVaR"] - extreme_values[0.0]["CVaR"])
        transform_sys = 1/(extreme_values[0.0]["System Expected"] - extreme_values[1.0]["System Expected"])
        new_budget = @constraint(model, c_budget_transformed, transform_sys*(model[:exp_sys_cost]) + transform_cvar*(model[:cvar_term]) <= budget)
        push!(budget_added, :c_budget_transformed)
        tranformed_dict = Dict("budget" => budget, "transform_cvar" => transform_cvar, "transform_sys" => transform_sys)
        existing_budgets[budget_type] = tranformed_dict
        constraint_dict[budget_type] = new_budget
    else
        error("Unsupported budget type. Please choose from 'Investment', 'Expected', 'CVaR', 'Transformed','System_Expected', or 'System_CVaR'.")
    end

    log_detail("Added budget constraint: $(model[budget_added[1]])")
    return existing_budgets, constraint_dict
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
    log_detail("Updated budget constraint ($budget_type) RHS to: $new_budget")
    return existing_budgets
end
