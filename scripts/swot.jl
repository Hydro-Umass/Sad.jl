### Run SAD algorithm with SWOT data

using Sad
using DelimitedFiles
using Distributions
using LinearAlgebra
using NCDatasets

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
= `swordfile`: SWORD NetCDF file

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
    write_output(reachid, valid, outdir, A0, n, Qa, Qu)

Write SAD output to NetCDF.

"""
function write_output(reachid, valid, outdir, A0, n, Qa, Qu, W, time_str)
    outfile = joinpath(outdir, "$(reachid)_sad.nc")
    out = Dataset(outfile, "c")
    out.attrib["valid"] = valid   # FIXME Determine what is considered valid in the context of a SAD run
    defDim(out, "nx", size(W, 1))
    defDim(out, "nt", size(W, 2))
    ridv = defVar(out, "reach_id", Int64, (), fillvalue = FILL)
    ridv[:] = reachid
    A0v = defVar(out, "A0", Float64, (), fillvalue = FILL)
    A0v[:] = A0
    nv = defVar(out, "n", Float64, (), fillvalue = FILL)
    nv[:] = n
    Qav = defVar(out, "Qa", Float64, ("nt",), fillvalue = FILL)
    Qav[:] = replace!(Qa, NaN=>FILL)
    Quv = defVar(out, "Q_u", Float64, ("nt",), fillvalue = FILL)
    Quv[:] = Qu
    time_str_var = defVar(out,"time_str", String, ("nt",), fillvalue = "no_data")
    time_str_var[:] = time_str
    close(out)
end

"""
    main()

Main driver routine.

"""
function main(reachid, swordfile, sosfile, swotfile, outdir)

    nids, x = river_info(reachid, swordfile)
    H, W, S, dA, Hr, Wr, Sr, time_str = read_swot_obs(swotfile, nids)

    reach = Sad.preprocess(x, H, W, S)

    A0 = missing
    n = missing
    Qa = Matrix{Sad.FloatM}(missing, 1, size(W, 2))
    Qu = Matrix{Sad.FloatM}(missing, 1, size(W, 2))
    if all(ismissing, H) || all(ismissing, W) || all(ismissing, S)
        println("$(reachid): INVALID")
        write_output(reachid, 0, outdir, A0, n, Qa, Qu, W, time_str)
    else
        p = Sad.priors(sosfile, reach.hmin, reachid)
        if ismissing(p)
            println("$(reachid): INVALID, missing mean discharge")
            write_output(reachid, 0, outdir, A0, n, Qa, Qu, W, time_str)
        else
            try
                res = Sad.infer(p, reach)
                A0  = Sad.compute_A0(reach, res.reach_ensemble)
                n   = mean(res.reach_ensemble[1, :])
                Qa[1, :]  = [isnan(q) ? missing : q for q in res.Q_post]
                Qu[1, :] = [isnothing(res.A_post[t]) ? missing : std(res.A_post[t][1, :]) for t in 1:reach.nt]
                println("$(reachid): VALID")
                write_output(reachid, 1, outdir, A0, n, Qa, Qu, W, time_str)
            catch
                println("$(reachid): INVALID")
                write_output(reachid, 0, outdir, A0, n, Qa, Qu, W, time_str)
            end
        end
    end
end
