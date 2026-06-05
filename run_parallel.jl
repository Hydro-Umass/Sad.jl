#!/usr/bin/env julia
#
# Run SAD v2 on all SWOT files in a directory in parallel.
#
# Usage:
#   julia -p 8 --project=. run_parallel.jl
#   julia -p 8 --project=. run_parallel.jl /data/swot
#   julia -p 8 --project=. run_parallel.jl /data/swot /output/v2
#   julia -p 8 --project=. run_parallel.jl --sword-dir /sword --sos-dir /sos
#   julia -p 8 --project=. run_parallel.jl --no-width   # disable width likelihood
#
# The -p N flag spawns N worker processes. Combine with -t N for threads:
#   julia -p 4 -t 4 --project=. run_parallel.jl
#
# All kwargs from main() are supported via flags:
#   --λ-smooth 0.1   --ν 5.0   --σ-obs NaN   --iterations 500   --g-tol 1e-6
#   --sword-dir DIR  --sos-dir DIR  --sword-version v17  --sos-version v17
#   --dry-run   --no-width   --dingman   --outlier-thresh 1.0
#

using Distributed
using Dates

# ── Simple argument parsing (no external deps) ─────────────────────────────

struct Config
    swot_dir::String
    sword_dir::String
    sos_dir::String
    out_dir::String
    sword_version::String
    sos_version::String
    λ_smooth::Float64
    ν::Float64
    σ_obs::Float64
    iterations::Int
    g_tol::Float64
    dry_run::Bool
    use_width::Bool
    rectangular::Bool
    outlier_thresh::Float64
end

function parse_args(args)
    defaults = Config(
        "../data/swot",          # swot_dir
        "../data/sword",         # sword_dir
        "../data/sos",           # sos_dir
        "../output",             # out_dir
        "v17",                   # sword_version
        "v17",                   # sos_version
        0.1,                     # λ_smooth
        5.0,                     # ν
        NaN,                     # σ_obs
        500,                     # iterations
        1e-6,                    # g_tol
        false,                   # dry_run
        false,                   # use_width
        true,                    # rectangular (default: rectangular channel)
        1.0                      # outlier_thresh [m]
    )

    positional = String[]
    i = 1
    kw = Dict{String, String}()
    dry_run = false
    no_width = false
    no_rect = false
    while i <= length(args)
        a = args[i]
        if a == "--dry-run"
            dry_run = true
            i += 1
        elseif a == "--no-width"
            no_width = true
            i += 1
        elseif a == "--dingman"
            no_rect = true
            i += 1
        elseif startswith(a, "--")
            key = a[3:end]  # strip leading "--"
            i += 1
            i > length(args) && error("Option $a requires a value")
            kw[key] = args[i]
            i += 1
        else
            push!(positional, a)
            i += 1
        end
    end

    swot_dir = get(kw, "swot-dir", get(positional, 1, defaults.swot_dir))
    out_dir  = get(kw, "out-dir", get(positional, 2, defaults.out_dir))

    use_width = no_width ? false : parse(Bool, get(kw, "use-width", string(defaults.use_width)))
    rectangular = no_rect ? false : parse(Bool, get(kw, "rectangular", string(defaults.rectangular)))

    Config(
        swot_dir,
        get(kw, "sword-dir", defaults.sword_dir),
        get(kw, "sos-dir", defaults.sos_dir),
        out_dir,
        get(kw, "sword-version", defaults.sword_version),
        get(kw, "sos-version", defaults.sos_version),
        parse(Float64, get(kw, "λ-smooth", string(defaults.λ_smooth))),
        parse(Float64, get(kw, "ν", string(defaults.ν))),
        parse(Float64, get(kw, "σ-obs", string(defaults.σ_obs))),
        parse(Int, get(kw, "iterations", string(defaults.iterations))),
        parse(Float64, get(kw, "g-tol", string(defaults.g_tol))),
        dry_run,
        use_width,
        rectangular,
        parse(Float64, get(kw, "outlier-thresh", string(defaults.outlier_thresh)))
    )
end

# Helper: indexed access with fallback
get(v::Vector, idx::Int, default) = idx <= length(v) ? v[idx] : default
get(d::Dict{String,String}, key::String, default) = haskey(d, key) ? d[key] : default

# ── Continent lookup (same convention as run_all.jl) ───────────────────────

function continent(id::Int)
    id10 = floor(Int, id / 1e10)
    return Dict(1 => "af", 2 => "eu", 3 => "as", 4 => "as", 5 => "oc",
                6 => "sa", 7 => "na", 8 => "na", 9 => "na")[id10]
end

# ── Build the task list from SWOT directory ─────────────────────────────────

function build_tasks(cfg::Config)
    swot_files = sort(filter(f -> endswith(f, "_SWOT.nc"),
                             readdir(cfg.swot_dir; join=true)))

    if isempty(swot_files)
        error("No *_SWOT.nc files found in '$(cfg.swot_dir)'")
    end

    tasks = NamedTuple[]
    for swotfile in swot_files
        reachid_str = replace(basename(swotfile), "_SWOT.nc" => "")
        reachid = parse(Int, reachid_str)
        cont = continent(reachid)

        swordfile = if isempty(cfg.sword_dir)
            "../data/sword/$(cont)_sword_$(cfg.sword_version).nc"
        else
            joinpath(cfg.sword_dir, "$(cont)_sword_$(cfg.sword_version).nc")
        end

        sosfile = if isempty(cfg.sos_dir)
            "../data/sos/$(cont)_sword_$(cfg.sos_version)_SOS_priors.nc"
        else
            joinpath(cfg.sos_dir, "$(cont)_sword_$(cfg.sos_version)_SOS_priors.nc")
        end

        push!(tasks, (reachid=reachid, swordfile=swordfile,
                      sosfile=sosfile, swotfile=swotfile))
    end

    return tasks
end

# ── Worker function (defined on all processes) ────────────────────────────

@everywhere function _run_sad_reach(reachid, swordfile, sosfile, swotfile, outdir;
                                    σ_obs, ν, λ_smooth, iterations, g_tol, use_width,
                                    rectangular, outlier_thresh)
    try
        main(reachid, swordfile, sosfile, swotfile, outdir;
             σ_obs=σ_obs, ν=ν, λ_smooth=λ_smooth,
             iterations=iterations, g_tol=g_tol,
             use_width=use_width, rectangular=rectangular,
             outlier_thresh=outlier_thresh)
        return (reachid=reachid, status=:ok)
    catch e
        @warn "$(reachid): unhandled error" exception=(e, catch_backtrace())
        return (reachid=reachid, status=:error, exception=string(e))
    end
end

# ── Main ───────────────────────────────────────────────────────────────────

function main_entry()
    cfg = parse_args(ARGS)
    tasks = build_tasks(cfg)

    println("=== SAD v2 Parallel Runner ===")
    println("SWOT directory:   $(abspath(cfg.swot_dir))")
    println("SWORD directory:   $(isempty(cfg.sword_dir) ? "../data/sword (default)" : abspath(cfg.sword_dir))")
    println("SoS directory:    $(isempty(cfg.sos_dir) ? "../data/sos (default)" : abspath(cfg.sos_dir))")
    println("Output directory:  $(abspath(cfg.out_dir))")
    println("Number of reaches: $(length(tasks))")
    println("Workers available:  $(nworkers())")
    println("Parameters: λ_smooth=$(cfg.λ_smooth), ν=$(cfg.ν), " *
            "σ_obs=$(cfg.σ_obs), iterations=$(cfg.iterations), g_tol=$(cfg.g_tol), " *
            "use_width=$(cfg.use_width), rectangular=$(cfg.rectangular), " *
            "outlier_thresh=$(cfg.outlier_thresh)m")

    if cfg.dry_run        println("\nDry run — would process:")
        for t in tasks
            cont = continent(t.reachid)
            println("  $(t.reachid) ($(cont)) → $(basename(t.swotfile))")
        end
        return
    end

    mkpath(cfg.out_dir)

    # Load swot.jl on all workers (which brings in Sad, NCDatasets, etc.)
    swot_script = joinpath(@__DIR__, "scripts", "swot.jl")
    @everywhere include($swot_script)

    # Extract plain values so the pmap closure only serializes
    # standard types — NOT the Config struct (workers don't have it).
    out_dir      = cfg.out_dir
    kw_σ_obs    = cfg.σ_obs
    kw_ν        = cfg.ν
    kw_λ_smooth = cfg.λ_smooth
    kw_iters    = cfg.iterations
    kw_g_tol    = cfg.g_tol
    kw_use_width = cfg.use_width
    kw_rectangular = cfg.rectangular
    kw_outlier_thresh = cfg.outlier_thresh

    println("\nProcessing $(length(tasks)) reaches on $(nworkers()) workers...")
    t_start = now()

    # pmap distributes one task per worker with automatic load balancing
    results = pmap(tasks) do task
        _run_sad_reach(task.reachid, task.swordfile, task.sosfile,
                       task.swotfile, out_dir;
                       σ_obs=kw_σ_obs, ν=kw_ν, λ_smooth=kw_λ_smooth,
                       iterations=kw_iters, g_tol=kw_g_tol,
                       use_width=kw_use_width, rectangular=kw_rectangular,
                       outlier_thresh=kw_outlier_thresh)
    end

    elapsed = canonicalize(now() - t_start)

    # Summary
    n_ok  = count(r -> r.status == :ok, results)
    n_err = count(r -> r.status == :error, results)

    println("\n=== Done in $(elapsed) ===")
    println("Succeeded: $n_ok / $(length(tasks))")
    if n_err > 0
        println("Failed:    $n_err / $(length(tasks))")
        for r in results
            if r.status == :error
                println("  ✗ reach $(r.reachid): $(r.exception)")
            end
        end
    end
end

main_entry()