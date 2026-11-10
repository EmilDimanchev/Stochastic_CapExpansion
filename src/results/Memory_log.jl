using Distributed, Printf

#==========================================================================================

Memory checkpoints for long cluster runs

log_memory_status(label) writes one line that starts with "MEM|" and then flushes stdout and
stderr. When Julia's output goes to a file (as under Slurm) both streams are block-buffered, so
without the flush the last part of the log is lost when the job is killed for running out of
memory, which is exactly the part that says where it was.

The line holds:
  - main process: current resident memory (RSS), its peak so far, and the Julia GC heap still in use (live)
  - workers: total and largest RSS and live heap, with the pid of the largest, and the largest peak RSS
    so far (a stage's peak only shows up here, since the checkpoints run between stages)
  - free memory on the machine (Sys.free_memory; under Slurm this is the node, not the job's limit)
  - main process GC time: cumulative, and its share of the wall time since the previous checkpoint

RSS is read from /proc/self/statm (Linux). Elsewhere (macOS) it falls back to Sys.maxrss(),
the peak rather than the current value, and the line says "peak".

Worker figures need one remote call per worker, so only call this between stages, when the
workers are idle: a worker in the middle of a long solve answers only when it's done.

==========================================================================================#

const MEMORY_LOG_LAST = Ref((wall_ns = UInt64(0), gc_ns = 0))

# (rss bytes, rss_is_peak, GC live bytes) of the calling process
function process_memory()
    rss, peak = Sys.maxrss(), true
    try
        fields = split(read("/proc/self/statm", String))
        rss = parse(Int, fields[2]) * ccall(:getpagesize, Cint, ())
        peak = false
    catch
    end
    return (rss = Int(rss), rss_is_peak = peak, peak_rss = Int(Sys.maxrss()), live = Int(Base.gc_live_bytes()))
end

gigabytes(bytes::Real) = bytes / 1024^3

# Memory of every worker, gathered in parallel. Empty when there are no workers.
function worker_memory()
    pids = nprocs() > 1 ? workers() : Int[]
    mems = Vector{Any}(undef, length(pids))
    @sync for (i, pid) in enumerate(pids)
        @async mems[i] = remotecall_fetch(process_memory, pid)
    end
    return pids, mems
end

function log_memory_status(label::AbstractString; include_workers::Bool = true)
    main = process_memory()
    pids, mems = include_workers ? worker_memory() : (Int[], Any[])

    now_ns, gc_ns = time_ns(), Base.gc_num().total_time
    prev = MEMORY_LOG_LAST[]
    gc_share = prev.wall_ns == 0 ? NaN : (gc_ns - prev.gc_ns) / max(now_ns - prev.wall_ns, 1)
    since_s = prev.wall_ns == 0 ? NaN : (now_ns - prev.wall_ns) / 1e9
    MEMORY_LOG_LAST[] = (wall_ns = now_ns, gc_ns = gc_ns)

    peak_tag(m) = m.rss_is_peak ? " peak" : ""
    main_rss = main.rss_is_peak ? @sprintf("rss peak %.2f GB", gigabytes(main.rss)) :
                                  @sprintf("rss %.2f GB (peak %.2f GB)", gigabytes(main.rss), gigabytes(main.peak_rss))
    parts = [@sprintf("MEM|%s|main %s, live %.2f GB", label, main_rss, gigabytes(main.live))]
    w_rss_total = w_rss_max = w_live_total = w_live_max = 0
    w_rss_pid = 0
    if !isempty(mems)
        rss = [m.rss for m in mems]
        live = [m.live for m in mems]
        w_rss_total, w_live_total = sum(rss), sum(live)
        i = argmax(rss)
        w_rss_max, w_rss_pid, w_live_max = rss[i], pids[i], maximum(live)
        w_peak_max = maximum(m.peak_rss for m in mems)
        push!(parts, @sprintf("%d workers: rss%s total %.2f GB, max %.2f GB (pid %d), largest peak so far %.2f GB; live total %.2f GB, max %.2f GB",
                              length(pids), peak_tag(mems[i]), gigabytes(w_rss_total), gigabytes(w_rss_max), w_rss_pid,
                              gigabytes(w_peak_max), gigabytes(w_live_total), gigabytes(w_live_max)))
    end
    push!(parts, @sprintf("free %.1f GB", gigabytes(Sys.free_memory())))
    push!(parts, isnan(gc_share) ? @sprintf("main GC %.1f s total", gc_ns / 1e9) :
                                   @sprintf("main GC %.1f s total, %.0f%% of the last %.0f s", gc_ns / 1e9, 100 * gc_share, since_s))
    @info(join(parts, "|"))
    flush(stdout); flush(stderr)
    return (main = main, workers = mems, pids = pids, worker_rss_total = w_rss_total, worker_rss_max = w_rss_max,
            worker_rss_max_pid = w_rss_pid, worker_live_total = w_live_total, worker_live_max = w_live_max, gc_share = gc_share)
end
