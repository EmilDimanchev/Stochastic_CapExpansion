using Printf, Test

"""
Peak memory and time of nearest_distance (src/model/sampling/gridfill_sampler.jl), full matrix vs
row chunks, at the sizes of the 2-zone Della GridFill run: m valid blends in one cell (at most
6 × 8000 = 48k) against n solved plans (about 555 when GridFill starts, about 2,555 after 2000
samples), d = 14 capacities.

Each case runs in a fresh Julia process and reports how much its peak resident memory (Sys.maxrss)
grew during the one call. @allocated can't be used: it adds up every block, so both versions
allocate about the same total. Distances are compared in the parent process.

Usage, from a REPL that already ran include("experiments/MGA_Tests.jl"):
  test_nearest_distance_memory()
The largest full-matrix case needs about 3 GB.
"""

isdefined(Main, :run_stochastic_exploration_risk_pareto) || include(joinpath(@__DIR__, "..", "experiments", "MGA_Tests.jl"))

const ND_REPO = abspath(joinpath(@__DIR__, ".."))

# Peak RSS growth (bytes) and time of one call, in a fresh process. variant: "full" or "chunked".
function nearest_distance_peak(variant::String, m::Int, n::Int; d::Int = 14, chunk::Int = 1024)
    script = """
    using LinearAlgebra, Random, JuMP
    include(joinpath($(repr(ND_REPO)), "src", "model", "sampling", "gridfill_bounds.jl"))
    include(joinpath($(repr(ND_REPO)), "src", "model", "sampling", "gridfill_sampler.jl"))
    full_nearest(Y, R, sc) = vec(minimum(pairwise_distances(Y ./ sc', R ./ sc'); dims = 2))
    f = $(repr(variant)) == "full" ? full_nearest : (Y, R, sc) -> nearest_distance(Y, R, sc; chunk = $chunk)
    rng = Xoshiro(3)
    sc = 1 .+ 100 .* rand(rng, $d)
    Y, R = 1000 .* rand(rng, $m, $d), 1000 .* rand(rng, $n, $d)
    f(Y[1:10, :], R[1:5, :], sc)          # compile
    GC.gc(); before = Sys.maxrss()
    t = @elapsed f(Y, R, sc)
    println(Sys.maxrss() - before, " ", t)
    """
    out = read(`$(Base.julia_cmd()) --project=$(ND_REPO) --startup-file=no -e $script`, String)
    grow, t = split(strip(split(strip(out), '\n')[end]))
    return parse(Int, grow), parse(Float64, t)
end

function nearest_distance_memory_table(; sizes = [(12_000, 555), (48_000, 555), (12_000, 2_555), (48_000, 2_555), (3_600, 2_555)], d = 14, chunk = 1024)
    full_nearest(Y, R, sc) = vec(minimum(pairwise_distances(Y ./ sc', R ./ sc'); dims = 2))
    println(@sprintf("%8s %7s | %13s %13s %8s | %8s %9s | %s", "m", "n", "full peak MiB", "chunked MiB", "ratio", "full s", "chunked s", "max |difference|"))
    rows = []
    for (m, n) in sizes
        p_full, t_full = nearest_distance_peak("full", m, n; d = d, chunk = chunk)
        p_chunk, t_chunk = nearest_distance_peak("chunked", m, n; d = d, chunk = chunk)
        rng = Xoshiro(3)
        sc = 1 .+ 100 .* rand(rng, d)
        Y, R = 1000 .* rand(rng, m, d), 1000 .* rand(rng, n, d)
        diff = maximum(abs.(nearest_distance(Y, R, sc; chunk = chunk) .- full_nearest(Y, R, sc)))
        println(@sprintf("%8d %7d | %13.1f %13.1f %7.1fx | %8.3f %9.3f | %.2e", m, n, p_full / 2^20, p_chunk / 2^20, p_full / max(p_chunk, 1), t_full, t_chunk, diff))
        push!(rows, (m = m, n = n, peak_full = p_full, peak_chunked = p_chunk, t_full = t_full, t_chunked = t_chunk, diff = diff))
        GC.gc()
    end
    return rows
end

@testset "nearest_distance peak memory, full vs chunked" begin
    rows = nearest_distance_memory_table()
    @test all(r.diff <= 1e-9 * 1000 for r in rows)
    # The chunked peak doesn't grow with m (a fixed allocator overhead of tens of MB sets its floor),
    # so the saving grows with the number of blends; the full version grows linearly in m.
    peak(m, n, field) = getfield(only(r for r in rows if r.m == m && r.n == n), field)
    for n in (555, 2_555)
        @test peak(48_000, n, :peak_chunked) <= 1.5 * peak(12_000, n, :peak_chunked)
        @test peak(48_000, n, :peak_full) >= 3 * peak(12_000, n, :peak_full)
        @test peak(48_000, n, :peak_full) >= 5 * peak(48_000, n, :peak_chunked)
    end
end
