using Logging, Printf, Dates, DataFrames

#==========================================================================================

Run log: what a run prints to stdout/stderr, and in what form

The "Log style" setting chooses between two loggers (configure_run_output!, called from
configure_parallel_workers!, which every experiment entry point calls first):

  "laptop" - everything, as before: every Benders iteration's progress lines, full JuMP
             expressions for objectives and budget constraints, and so on.
  "della"  - one line per message, stamped with the time since the run started
             ([+HH:MM:SS]) and flushed straight away, so an OOM-killed job keeps its log up to
             the last finished iteration (see the CLAUDE.md trap on lost .out tails). The
             per-iteration detail lines are hidden; ITER| lines carry what matters from them.

Detail-only messages go through log_detail, at the DETAIL level (between Debug and Info): the
laptop logger shows them labelled as Info, the della logger drops them. Package @debug
messages stay hidden in both styles.

Lines shared by both styles, written to be read by a person and parsed by a script
(parse_run_log reads ITER| and SOLVE| back into DataFrames):

  RUN| ...        run header: start time, style, workers, scenarios, run settings
  === ...         phase banners: anchors, each MGA vector, GridFill, final write, done
  ITER| label | it j | gap g% | LB a UB b | MP t s (sep +n in r rnd) | [reg t s |] SP t s | live n
                  one per Benders iteration, describing iterate j: the MP solve that produced it
                  (with separation rounds), the regularization solve if one ran, and the
                  sub-problem round that evaluated it. MGA solves show one gap and UB per budget
                  type ("gap T 0.05% CVaR 0.16% | UB T 51.12 CVaR 1.01e5"); their LB is the
                  budget RHS on the SOLVE| line. `live` is the MP's live cut count at that point.
  SOLVE| tag | label | [budget ... |] N it, gap g% | wall t s (MP, [reg,] SP, final) | cuts ... | costs
                  one per converged solve, after its ITER| lines. `final` is the work after
                  convergence (MP re-solve, full-payload sub-problems, eval SPs), which the ITER|
                  times don't include. Costs are $B per year.
  MEM| ...        memory checkpoints (Memory_log.jl)
  PHASE| ...      della only: written (and flushed) right before each step of the master problem that
                  calls into Gurobi (seed, solve, archive delete, trim, level-set solve), so after a crash
                  in native code (no Julia backtrace) the last line names the step it died in
  TRIMRUN| ...    cut-trim diagnostics, only when a Run log vector is passed

==========================================================================================#

const DETAIL = LogLevel(-500)
const RUN_LOG_STYLES = ("laptop", "della")

# A message only the laptop style shows. A function rather than a macro: src/ files are
# included alphabetically, so model/ code couldn't use a macro defined in results/.
log_detail(msg::AbstractString) = @logmsg DETAIL msg

# True when the current logger shows DETAIL messages (the laptop style)
detail_enabled() = Logging.min_enabled_level(current_logger()) <= DETAIL

# A PHASE| marker (see the header). Della logger only: it flushes every line it writes, and the
# laptop style already shows the detail lines around each of these steps. (CompactLogger is defined
# below; the check runs at call time.)
log_phase(msg::AbstractString) = current_logger() isa CompactLogger && @info("PHASE| " * msg)

# "Log style" from settings; settings files written before the key existed get "laptop"
function run_log_style(settings::AbstractDict)
    style = get(settings, "Log style", "laptop")
    style in RUN_LOG_STYLES || error("\"Log style\" must be one of $(RUN_LOG_STYLES), got $(repr(style)).")
    return style
end

# Laptop: the default console format, with DETAIL labelled as Info so lines look as they always have
laptop_metafmt(level::LogLevel, args...) = Logging.default_metafmt(level == DETAIL ? Logging.Info : level, args...)

#==========================================================================================
Della logger
==========================================================================================#

struct CompactLogger <: AbstractLogger
    stream::IO
    min_level::LogLevel
    start::Float64 # time() the elapsed stamps count from
end
CompactLogger(stream::IO = stderr; min_level::LogLevel = Logging.Info, start::Float64 = time()) = CompactLogger(stream, min_level, start)

Logging.min_enabled_level(logger::CompactLogger) = logger.min_level
Logging.shouldlog(::CompactLogger, args...) = true
Logging.catch_exceptions(::CompactLogger) = true

function elapsed_stamp(seconds::Real)
    s = max(floor(Int, seconds), 0)
    return @sprintf("[+%02d:%02d:%02d]", s ÷ 3600, (s % 3600) ÷ 60, s % 60)
end

# Plain Info messages print as one stamped line. Warnings, errors and anything with key-value
# pairs (e.g. exception = (e, backtrace)) are rendered by a ConsoleLogger, so nothing is lost,
# and their continuation lines are indented under the stamp.
function Logging.handle_message(logger::CompactLogger, level, message, _module, group, id, file, line; kwargs...)
    stamp = elapsed_stamp(time() - logger.start)
    text = if level < Logging.Warn && isempty(kwargs)
        string(message)
    else
        buf = IOBuffer()
        Logging.handle_message(ConsoleLogger(buf, Logging.BelowMinLevel), level, message, _module, group, id, file, line; kwargs...)
        rstrip(String(take!(buf)))
    end
    pad = " "^(length(stamp) + 1)
    out = IOBuffer()
    for (i, ln) in enumerate(split(text, '\n'))
        print(out, i == 1 ? stamp * " " : pad, ln, '\n')
    end
    write(logger.stream, take!(out))
    flush(logger.stream)
    flush(stdout) # println output (e.g. @time) is block-buffered too when it goes to a file
    return nothing
end

function run_logger(style::AbstractString, stream::IO = stderr)
    style == "della" && return CompactLogger(stream)
    return ConsoleLogger(stream, DETAIL; meta_formatter = laptop_metafmt)
end

# Sets the global logger for this run from "Log style". Returns the style.
function configure_run_output!(settings::AbstractDict; stream::IO = stderr)
    style = run_log_style(settings)
    global_logger(run_logger(style, stream))
    return style
end

#==========================================================================================
Line builders (pure functions, so they can be tested without a run)
==========================================================================================#

const BUDGET_ABBREVIATIONS = Dict("Transformed" => "T", "System_Expected" => "SysEV", "Expected" => "EV", "CVaR" => "CVaR",
                                  "System_Weighted_CVaR" => "SysWCVaR", "Investment" => "Inv")
budget_abbreviation(budget_type::AbstractString) = get(BUDGET_ABBREVIATIONS, budget_type, budget_type)

fmt_pct(fraction::Real) = @sprintf("%.3g%%", 100 * fraction)
fmt_val(x::Real) = @sprintf("%.7g", x)
fmt_secs(x::Real) = @sprintf("%.3f", x) # ms resolution: laptop MP solves take a few ms
fmt_count(n::Integer) = n < 1000 ? string(n) : @sprintf("%.1fk", n / 1000)

# `gap` and `UB` are numbers for a single-objective (anchor) solve, or vectors of
# budget type => value pairs for an MGA solve. `sep` is the separate! stats of the MP solve
# that produced this iterate; `time_reg` is nothing when no regularization solve ran.
function iter_line(label::AbstractString, j::Integer; gap, UB, LB = nothing, time_mp::Real, time_sp::Real, sep, live::Integer, time_reg = nothing)
    parts = ["ITER| $label", "it $j"]
    if gap isa Real
        push!(parts, "gap $(fmt_pct(gap))")
        push!(parts, LB === nothing ? "UB $(fmt_val(UB))" : "LB $(fmt_val(LB)) UB $(fmt_val(UB))")
    else
        push!(parts, "gap " * join(["$(budget_abbreviation(t)) $(fmt_pct(g))" for (t, g) in gap], " "))
        push!(parts, "UB " * join(["$(budget_abbreviation(t)) $(fmt_val(u))" for (t, u) in UB], " "))
    end
    push!(parts, "MP $(fmt_secs(time_mp)) s (sep +$(sep.restored) in $(sep.rounds) rnd)")
    time_reg === nothing || push!(parts, "reg $(fmt_secs(time_reg)) s")
    push!(parts, "SP $(fmt_secs(time_sp)) s")
    push!(parts, "live $live")
    return join(parts, " | ")
end

# The budget right-hand sides active for a solve, from the `budgets` Dict the drivers keep
function budget_summary(budgets::AbstractDict)
    entries = String[]
    for (t, v) in sort!(collect(budgets); by = first)
        rhs = v isa AbstractDict ? v["budget"] : v
        push!(entries, "$(budget_abbreviation(t))≤$(@sprintf("%.6g", rhs))")
    end
    return "budget " * join(entries, " ")
end

# One converged solve. `out` is the Dict benders_algorithm / mga_benders return; `weight` is the
# anchor's risk weight (adds the risk-adjusted cost); `pool_size` is live + archived cuts.
function solve_line(tag::AbstractString, label::AbstractString, out::AbstractDict; pool_size::Integer, budgets = nothing, weight = nothing)
    gaps = out["Gaps"]
    n_iter = length(out["Time SP hist"])
    times = ["MP $(@sprintf("%.1f", sum(out["Time MP hist"])))"]
    haskey(out, "Time Reg hist") && push!(times, "reg $(@sprintf("%.1f", sum(out["Time Reg hist"])))")
    push!(times, "SP $(@sprintf("%.1f", sum(out["Time SP hist"])))", "final $(@sprintf("%.1f", out["Time final"]))")
    restored = sum((s.restored for s in out["Sep hist"]); init = 0)

    inv, ev, cvar = out["MP"]["Inv_cost"] / 1e9, out["Expected Value"] / 1e9, out["CVaR"] / 1e9
    costs = @sprintf("inv %.4g EV %.4g CVaR %.4g", inv, ev, cvar)
    weight === nothing || (costs *= @sprintf(" risk-adj %.4g", inv + weight * ev + (1 - weight) * cvar))

    parts = ["SOLVE| $tag", label]
    budgets === nothing || isempty(budgets) || push!(parts, budget_summary(budgets))
    push!(parts, "$n_iter it, gap $(fmt_pct(maximum(gaps[end])))")
    push!(parts, "wall $(@sprintf("%.1f", out["Elapsed"])) s ($(join(times, ", ")))")
    push!(parts, "cuts $(fmt_count(out["Live cuts at start"])) live at start, +$(fmt_count(restored)) restored, pool $(fmt_count(pool_size))")
    push!(parts, costs * " (\$B/yr)")
    return join(parts, " | ")
end

# First RUN| line: what every run shares (the drivers add a second line with their own settings)
function run_header(settings::AbstractDict, inputs::AbstractDict)
    S, F, K = inputs["Number of demand scenarios"], inputs["Number of gas price scenarios"], inputs["Number of weather scenarios"]
    B = get(settings, "Subproblems per scenario", 1)
    workers = settings["Parallel flag"] && nprocs() > 1 ? "$(nworkers()) workers" : "serial"
    return "RUN| started $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS")) | log style $(run_log_style(settings)) | $workers | " *
           "$(S * F * K) scenarios ($S demand × $F gas × $K weather) × $B block$(B == 1 ? "" : "s") | " *
           "tol $(fmt_pct(settings["Convergence tolerance"])), max $(max_benders_iterations(settings)) it | solver $(settings["Solver"]) | costs in \$B/yr"
end

function fmt_duration(seconds::Real)
    s = round(Int, seconds)
    return s < 3600 ? @sprintf("%dm%02ds", s ÷ 60, s % 60) : @sprintf("%dh%02dm", s ÷ 3600, (s % 3600) ÷ 60)
end

# Conflicting constraints from an infeasible MP. Laptop: the full constraints, as before. Della:
# names only (bounds, which have no name, as their short expression), at most `max_shown`.
function show_conflicts(cons::AbstractVector; max_shown::Integer = 25)
    if detail_enabled()
        display(cons)
        return nothing
    end
    labels = [isempty(name(c)) ? string(c) : name(c) for c in cons]
    shown = join(first(labels, max_shown), ", ")
    length(labels) > max_shown && (shown *= ", … ($(length(labels) - max_shown) more)")
    @info("Conflicting constraints ($(length(cons))): $shown")
    return nothing
end

#==========================================================================================
Parser: ITER| and SOLVE| lines of a .out file back into tables
==========================================================================================#

const STAMP_RE = r"\[\+(\d+):(\d\d):(\d\d)\]"

function _stamp_seconds(line::AbstractString)
    m = match(STAMP_RE, line)
    return m === nothing ? missing : 3600 * parse(Int, m[1]) + 60 * parse(Int, m[2]) + parse(Int, m[3])
end

_num(s::AbstractString) = parse(Float64, s)
_count(s::AbstractString) = endswith(s, "k") ? round(Int, 1000 * _num(chop(s))) : parse(Int, s)
_max_pct(s::AbstractString) = maximum(_num(m[1]) for m in eachmatch(r"(-?[\d.]+(?:e[-+]?\d+)?)%", s))

# Returns (iters = DataFrame, solves = DataFrame). Works on either log style (laptop lines have no
# elapsed stamp, so `elapsed_s` is missing there). Gaps are in percent, times in seconds.
function parse_run_log(path::AbstractString)
    iters = DataFrame(elapsed_s = Union{Missing,Int}[], label = String[], iteration = Int[], gap_pct = Float64[], mp_s = Float64[],
                      sep_restored = Int[], sep_rounds = Int[], reg_s = Union{Missing,Float64}[], sp_s = Float64[], live = Int[])
    solves = DataFrame(elapsed_s = Union{Missing,Int}[], tag = String[], label = String[], iterations = Int[], gap_pct = Float64[],
                       wall_s = Float64[], mp_s = Float64[], reg_s = Union{Missing,Float64}[], sp_s = Float64[], final_s = Float64[])
    for line in eachline(path)
        if (i = findfirst("ITER| ", line)) !== nothing
            parts = split(line[last(i)+1:end], " | ")
            field = prefix -> something(findfirst(p -> startswith(p, prefix), parts), 0)
            mp = match(r"MP ([\d.]+) s \(sep \+(\d+) in (\d+) rnd\)", parts[field("MP ")])
            reg = field("reg ") == 0 ? missing : _num(match(r"reg ([\d.]+) s", parts[field("reg ")])[1])
            push!(iters, (_stamp_seconds(line), String(parts[1]), parse(Int, split(parts[2])[2]), _max_pct(parts[3]),
                          _num(mp[1]), parse(Int, mp[2]), parse(Int, mp[3]), reg,
                          _num(match(r"SP ([\d.]+) s", parts[field("SP ")])[1]), parse(Int, split(parts[field("live ")])[2])))
        elseif (i = findfirst("SOLVE| ", line)) !== nothing
            parts = split(line[last(i)+1:end], " | ")
            it_part = parts[findfirst(p -> occursin(r"^\d+ it, gap", p), parts)]
            wall_part = parts[findfirst(p -> startswith(p, "wall "), parts)]
            t = key -> (m = match(Regex("\\b$key ([\\d.]+)"), wall_part)) === nothing ? missing : _num(m[1])
            push!(solves, (_stamp_seconds(line), String(parts[1]), String(parts[2]), parse(Int, split(it_part)[1]), _max_pct(it_part),
                           _num(match(r"wall ([\d.]+) s", wall_part)[1]), t("MP"), t("reg"), t("SP"), t("final")))
        end
    end
    return (iters = iters, solves = solves)
end
