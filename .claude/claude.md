### Running instructions for key functions in this repo:

# How to open julia environment/run functions in Stochastic_CapExploration:
1. change pwd to Stochastic_CapExploration
2. run <julia --project=.>
3. run <include("experiments/MGA_Tests.jl")>
4. run your desired function with test_index set to a name

# General instructions:
 - Please make tests for your functions, and make them runable with a similar command structure as above
 - Please only use the "laptop" variant of functions using inputs from Stochastic_CapExpansion/inputs/Inputs_30d_1000scen_7tech_2z
 - If a rapid test is needed, please feel free to disable one of the uncertainty sets for the model construction. Please leave one on though to maintain stochasticity.

# Traps
Recurring mistakes (Claude- or user-side) worth avoiding next time. Add to this list whenever a new one is caught.

- **No `timeout` shell command on this machine.** macOS doesn't ship GNU coreutils' `timeout`
  (would need `gtimeout` via `brew install coreutils`). Don't prefix Julia/other long-running
  commands with `timeout N ...` - it just fails with "command not found". Use the Bash tool's own
  `timeout` parameter instead.
- **`include("src/Stochastic_CapExpansion.jl")` fails standalone.** `src/model/benders/
  subproblems_economic_dispatch.jl` and `master_planning.jl` use `@fetchfrom`, `@sync`,
  `remotecall_wait`, `workers()`, etc. without importing `Distributed` themselves - they rely on the
  caller having done `using Distributed` first (exactly what `experiments/MGA_Tests.jl` does on its
  first line). Loading the module any other way (e.g. a quick REPL check) needs `using Distributed`
  run *before* the `include`, or it fails with `UndefVarError: @fetchfrom not defined`.
- **`set_optimizer(model, ...)` silently discards every previously-set attribute**, not just the one
  you meant to change - it swaps in a completely fresh backend with default settings. Any code path
  that reads "restore the optimizer" after a temporary solver swap (e.g. HiGHS -> Clarabel/Ipopt and
  back, see `master_planning.jl`'s regularization functions) needs to re-specify every attribute the
  base build set (tolerances, `run_crossover`, `Method`, etc.) via `optimizer_with_attributes`/
  `set_optimizer_attribute`, not just the one attribute that motivated the swap - a bare
  `HiGHS.Optimizer()` with no attributes was found silently reverting tolerances that had been set
  elsewhere. Grep for `set_optimizer(` whenever reasoning about what solver state persists across an
  iteration.
- **Multiple versions of a package can be installed side-by-side in `~/.julia/packages/<Pkg>/`**
  (different content-hash directories). Never assume the first one found is "the" version in use -
  cross-reference the `git-tree-sha1` in `Manifest.toml` (or just grep all installed versions for the
  API in question, since behavior is usually consistent across nearby versions) before citing solver
  API behavior as fact.
- **A commented-out "simpler fix" sitting next to a heavier workaround is not necessarily an
  unexplored idea.** It can be a deliberately abandoned dead end from prior empirical testing that
  didn't hold up. Don't propose reviving it from static code reading alone - ask the user first,
  especially for numerically-motivated workarounds that carry an explanatory comment (e.g. the
  HiGHS-\>Clarabel swap in `level_set_regularization`, which exists because the regularization
  problem needs a true interior-point solve to behave well, not because nobody tried crossover-on).
- **`include(...)` inside a script resolves relative paths against that script's own directory, not
  the process's cwd.** A one-off script living outside the repo (e.g. in a scratchpad) that does
  `include("experiments/MGA_Tests.jl")` will fail with "No such file or directory" even when run with
  `julia --project=. /path/to/script.jl` from the repo root - `--project=.` only sets the active
  project, it doesn't affect where relative `include` paths resolve. Use an absolute path (e.g.
  `joinpath(REPO_ROOT, "experiments", "MGA_Tests.jl")`) in any script that doesn't live in the repo
  root itself.
- **Don't trust a single post-change run's internal trend as proof of a performance win.** A solve
  time that drops sharply across iterations in the *new* code can be present in the *unmodified* code
  too (e.g. because per-call overhead dominates at small problem sizes, unrelated to the change made).
  Always stash the change and re-run the identical scenario on the original code for a true baseline
  before reporting a speedup - not just eyeballing the trend within one run.
- **The "laptop" input config (`Inputs_30d_1000scen_7tech_2z`) is too small to measure solver-level
  timing differences** (warm-start settings, `Method`/`Crossover` choices, etc.) - individual
  `optimize!` calls complete in single-digit-to-low-double-digit milliseconds, likely dominated by
  fixed Julia/JuMP/MOI call overhead rather than actual solver work. It's the right config for
  *correctness* checks (objective/gap/capacity-mix should match exactly across such changes) but not
  for demonstrating a performance improvement - that needs a larger-scale input set, which doesn't
  exist in this repo yet.
