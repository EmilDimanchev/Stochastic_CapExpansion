#-------------------------------------------------------------------------------------------
# Settings dicts for run_stochastic_exploration_risk_pareto, one per subsection of the run.
# Callers pass only the keys they want to change; everything else takes the default below.
# Unknown keys error out, so a typo can't silently fall back to a default.
#-------------------------------------------------------------------------------------------

# Budget constraint used during capacity exploration (and its optional tightening schedule)
default_budget_settings() = Dict{String,Any}(
    "Budget type" => "Transformed",       # "Transformed" or "Box"
    "Budget offset" => 1.0,               # Transformed: added to the loosest budget
    "Budget multiplier" => 1.10,          # Box: scales the EV-optimal System Expected cost
    "Tighten budget" => false,            # Transformed only: step the budget down per MGA vector
    "Tighten step size" => 0.1,
    "Floor offset" => 0.01,               # tightest budget = min transformed anchor cost + this
)

# Extra caps on one metric (System Expected / CVaR) at selected tightening levels
default_cap_settings() = Dict{String,Any}(
    "Cap wings" => true,                  # only active when tightening a Transformed budget
    "Number of caps" => 2,                # max caps per wing per level
    "Cap levels" => [1],                  # tightening levels at which caps are run
)

# Mapping of the feasible space (records every Benders iterate, then samples the interior)
default_mapping_settings() = Dict{String,Any}(
    "Mapping" => false,
    "Sampling" => false,
    "Interior samples" => 100,            # number of interior portfolios to evaluate (either method)
    "Sampling method" => "Delaunay",      # "Delaunay" (cost-space simplices) or "Boundary" (see below)
)

# Boundary interior sampler (src/model/boundary_sampling): picks blends of neighbouring final
# solves just inside the hull boundary, spread out by greedy max-min distance in a metric
# built from the Benders capacity gradients. Only used when "Sampling method" => "Boundary".
# Shared choices live elsewhere, as for the Delaunay sampler: the number of samples is
# mapping_settings["Interior samples"] and the RNG seed is settings["Seed"] (settings.yaml).
default_boundary_settings() = Dict{String,Any}(
    "Quotas" => ["outer_shell" => 0.50, "ev_wing" => 0.15, "cvar_wing" => 0.15, "interior" => 0.20], # base-region shares, in selection order
    "K neighbours" => 6,                  # neighbours a base vertex can blend with
    "Base weight range" => (0.5, 0.85),   # weight kept on the base vertex
    "Candidates per pick" => 30,          # candidate pool size per pick
    "Metric" => "active_subspace",        # "active_subspace" or "lengthscale_file"
    "Metric file" => nothing,             # length-scale CSV (tech, cap_max_MW, length_scale); only read for "lengthscale_file"
    "Energy share" => 0.999,              # keep the fewest eigen-directions reaching this share (0.95 kept r = 1 on run 13 and over-spread the picks)
    "r max" => 50,                        # ... but at most this many
    "Explicit M max dim" => 200,          # d above this: thin SVD instead of forming the d x d matrix
    "Waves" => 1,                         # split the interior samples into this many select/evaluate rounds
    "Refit metric each wave" => false,    # rebuild the subspace with the interior gradients included
    "Write gradients" => false,           # also write per-solve Gradients.csv files to disk
)

# Run bookkeeping: MGA vectors, output naming, evaluation SPs and diagnostics
default_run_settings() = Dict{String,Any}(
    "Vector set" => nothing,              # nothing -> generate from settings["Vector Type"]
    "Summary name" => "new_setup",
    "Eval SPs" => nothing,
    "Run log" => nothing,                 # Vector to collect per-run diagnostics (see log_run!)
)

function resolve_settings(defaults::Dict{String,Any}, overrides::AbstractDict, name::String)
    unknown = setdiff(keys(overrides), keys(defaults))
    isempty(unknown) || error("Unknown key(s) $(collect(unknown)) in $name. Valid keys: $(sort(collect(keys(defaults)))).")
    return merge(defaults, Dict{String,Any}(overrides))
end
