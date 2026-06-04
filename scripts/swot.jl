### Run SAD v2 algorithm with SWOT data

using Sad
using LinearAlgebra
using NCDatasets
using Statistics
using Optim

const FILL = -999999999999

"""
    read_swot_obs(ncfile, nids)

Load SWOT observations.

# Arguments

- `ncfile`: NetCDF file with SWOT observations
- `nids`: sorted node IDs from downstream to upstream

"""
function read_swot_obs(ncfile::String, nids::Vector{Int})
    Dataset(ncfile) do ds
        nodes = ds.group["node"]
        reaches = ds.group["reach"]
        S = permutedims(nodes["slope2"][:, :])
        H = permutedims(nodes["wse"][:, :])
        W = permutedims(nodes["width"][:, :])
        dA = reaches["d_x_area"][:]
        dA = convert(Vector{Sad.FloatM}, dA)
        Hr = convert(Vector{Sad.FloatM}, reaches["wse"][:])
        Wr = convert(Vector{Sad.FloatM}, reaches["width"][:])
        Sr = convert(Vector{Sad.FloatM}, reaches["slope2"][:])
        W[.!ismissing.(H) .&& isnan.(H)] .= missing
        S[.!ismissing.(H) .&& isnan.(H)] .= missing
        H[.!ismissing.(H) .&& isnan.(H)] .= missing
        nid = nodes["node_id"][:]
        dmap = Dict(nid[k] => k for k=1:length(nid))
        i = [dmap[k] for k in nids]
        time = reaches["time"][:]
        time_str = [string(t) for t in time]

        H[i, :], W[i, :], S[i, :], dA, Hr, Wr, Sr, time_str
    end
end

"""
    river_info(id, swordfile)

Retrieve information about river reach cross sections.

# Arguments

- `id`: reach ID
- `swordfile`: SWORD NetCDF file

"""
function river_info(id::Int, swordfile::String)
    Dataset(swordfile) do fd
        g = fd.group["nodes"]
        i = findall(g["reach_id"][:] .== id)
        nid = g["node_id"][i]
        x = g["dist_out"][i]
        k = findall(.!ismissing.(nid))
        x = x[k]
        nid = nid[k]
        # subtract the minimum `dist_out` and then sort from downstream to upstream
        x = x .- minimum(x)
        j = sortperm(x)
        close(fd)
        convert(Vector{Int}, nid[j]), x[j]
    end
end

"""
    write_output(reachid, valid, outdir, res, nt, W, time_str)

Write SAD v2 output to NetCDF.

Saves MAP parameter estimates (Q, n, r, z₀) and Laplace uncertainty (Q_std).

# Arguments

- `reachid`:  SWORD reach ID
- `valid`:    whether inference succeeded (1 = valid, 0 = invalid)
- `outdir`:   output directory
- `res`:      result NamedTuple from `infer` (or `nothing` if invalid)
- `nt`:       number of timesteps
- `W`:        width observations (nodes × timesteps)
- `time_str`: timestamp strings for each timestep
"""
function write_output(reachid, valid, outdir, res, nt, W, time_str)
    outfile = joinpath(outdir, "$(reachid)_sad.nc")
    out = Dataset(outfile, "c")
    out.attrib["sad_version"] = "v2"
    out.attrib["valid"] = valid
    defDim(out, "nx", size(W, 1))
    defDim(out, "nt", size(W, 2))
    ridv = defVar(out, "reach_id", Int64, (), fillvalue = FILL)
    ridv[:] = reachid

    if valid == 1 && !isnothing(res)
        # MAP channel parameters
        nv = defVar(out, "n", Float64, (), fillvalue = FILL)
        nv[:] = Float64(res.n_post)
        rv = defVar(out, "r", Float64, (), fillvalue = FILL)
        rv[:] = Float64(res.r_post)
        z0v = defVar(out, "z0", Float64, (), fillvalue = FILL)
        z0v[:] = Float64(res.z0_post)

        # Discharge: posterior mean and uncertainty
        Qa = defVar(out, "Qa", Float64, ("nt",), fillvalue = FILL)
        Qa[:] = replace(res.Q_post, NaN => FILL)
        Qu = defVar(out, "Q_u", Float64, ("nt",), fillvalue = FILL)
        Qu[:] = replace(res.Q_std, NaN => FILL)

        # Convergence / fallback info
        fallback = haskey(res, :fallback) ? res.fallback : false
        out.attrib["fallback"] = fallback ? 1 : 0
        if !isnothing(res.result)
            out.attrib["converged"] = string(res.result.stopped_by.g_converged)
            out.attrib["iterations"] = res.result.iterations
            out.attrib["nll_minimum"] = Float64(Optim.minimum(res.result))
        elseif fallback
            out.attrib["converged"] = "false_fallback"
            out.attrib["iterations"] = 0
            out.attrib["nll_minimum"] = Float64(NaN)
        end
    else
        nv = defVar(out, "n", Float64, (), fillvalue = FILL)
        nv[:] = FILL
        rv = defVar(out, "r", Float64, (), fillvalue = FILL)
        rv[:] = FILL
        z0v = defVar(out, "z0", Float64, (), fillvalue = FILL)
        z0v[:] = FILL
        Qa = defVar(out, "Qa", Float64, ("nt",), fillvalue = FILL)
        Qa[:] = fill(FILL, nt)
        Qu = defVar(out, "Q_u", Float64, ("nt",), fillvalue = FILL)
        Qu[:] = fill(FILL, nt)
    end

    time_str_var = defVar(out, "time_str", String, ("nt",), fillvalue = "no_data")
    time_str_var[:] = time_str
    close(out)
end

"""
    main(reachid, swordfile, sosfile, swotfile, outdir; kwargs...)

Main driver routine for SAD v2 discharge estimation.

# Arguments

- `reachid`:   SWORD reach ID
- `swordfile`: path to SWORD NetCDF database
- `sosfile`:   path to SoS NetCDF prior database
- `swotfile`:  path to SWOT observation NetCDF
- `outdir`:    output directory

# Keyword arguments

- `σ_obs`:       observation error scale [m] (default NaN: auto-estimate)
- `ν`:           Student-t d.f.; `Inf` = Gaussian (default 5.0)
- `λ_smooth`:    temporal smoothness weight (default 0.1)
- `iterations`:  max L-BFGS iterations (default 500)
- `g_tol`:       gradient convergence tolerance (default 1e-6)
- `use_width`:   use width observations in likelihood (default true)
"""
function main(reachid, swordfile, sosfile, swotfile, outdir;
              σ_obs::Float64     = NaN,
              ν::Float64         = 5.0,
              λ_smooth::Float64  = 0.1,
              iterations::Int    = 500,
              g_tol::Float64     = 1e-6,
              use_width::Bool    = true)

    nids, x = river_info(reachid, swordfile)
    H, W, S, dA, Hr, Wr, Sr, time_str = read_swot_obs(swotfile, nids)

    if all(ismissing, H) || all(ismissing, W) || all(ismissing, S)
        println("$(reachid): INVALID — no observations")
        write_output(reachid, 0, outdir, nothing, size(W, 2), W, time_str)
        return
    end

    reach = Sad.preprocess(x, H, W, S)

    p = Sad.priors(sosfile, reach.hmin, reachid; S0=mean(reach.S0.(reach.x)))
    if ismissing(p)
        println("$(reachid): INVALID — missing mean discharge prior")
        write_output(reachid, 0, outdir, nothing, reach.nt, W, time_str)
        return
    end

    try
        res = Sad.infer(p, reach;
                            time_str  = time_str,
                            σ_obs     = σ_obs,
                            ν         = ν,
                            λ_smooth  = λ_smooth,
                            iterations = iterations,
                            g_tol     = g_tol,
                            use_width = use_width)

        # Check if inference produced valid results
        # (data-sparse reaches return NaN with fallback=true)
        has_result = !isnothing(res.result)
        is_fallback = haskey(res, :fallback) ? res.fallback : false
        is_valid = !isnan(res.n_post) && !is_fallback

        if is_valid
            n_valid = length(findall(reach.valid))
            println("$(reachid): VALID (n=$(round(res.n_post, digits=4)), " *
                    "r=$(round(res.r_post, digits=2)), " *
                    "z0=$(round(res.z0_post, digits=2)), " *
                    "$(n_valid) valid ts, " *
                    "converged=$(has_result ? res.result.stopped_by.g_converged : "N/A")")
            write_output(reachid, 1, outdir, res, reach.nt, W, time_str)
        else
            println("$(reachid): INVALID — insufficient data or inference failed")
            write_output(reachid, 0, outdir, nothing, reach.nt, W, time_str)
        end

    catch e
        @warn "$(reachid): v2 inference failed" exception=(e, catch_backtrace())
        println("$(reachid): INVALID — v2 inference error")
        write_output(reachid, 0, outdir, nothing, reach.nt, W, time_str)
    end
end