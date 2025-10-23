### Run SAD algorithm with SWOT data

using Sad
using DelimitedFiles
using Distributions
using LinearAlgebra
using NCDatasets
using JSON

const FILL = -999999999999

"""
    get_reach_files(indir, reachjson)

Get reach file names.

"""
function get_reach_files(indir, reachjson)
    line = try parse(Int64, ENV["AWS_BATCH_JOB_ARRAY_INDEX"]) + 1 catch KeyError 1 end
    open(joinpath(indir, reachjson)) do jf
        data = read(jf, String)
        reachlist = JSON.parse(data)[line]
        reachlist["reach_id"], joinpath(indir, "swot", reachlist["swot"]), joinpath(indir, "sos", reachlist["sos"]), joinpath(indir, "sword", reachlist["sword"])
    end
end

"""
    read_swot_obs(ncfile, nids)

Load SWOT observations.

# Arguments

- `ncfile`: NetCDF file with SWOT observations
- `nids`: sorted node IDs from downstream to upstream

"""
function read_swot_obs(ncfile::String, nids::Vector{Int})
    Dataset(ncfile) do ds
        nodes = NCDatasets.group(ds, "node")
        reaches = NCDatasets.group(ds, "reach")
        S = permutedims(nodes["slope2"][:])
        H = permutedims(nodes["wse"][:])
        W = permutedims(nodes["width"][:])
        dA = reaches["d_x_area"][:]
        dA = convert(Vector{Sad.FloatM}, dA)
        Hr = convert(Vector{Sad.FloatM}, reaches["wse"][:])
        Wr = convert(Vector{Sad.FloatM}, reaches["width"][:])
        Sr = convert(Vector{Sad.FloatM}, reaches["slope2"][:])
        W[.!ismissing.(H) .&& isnan.(H)] .= missing
        S[.!ismissing.(H) .&& isnan.(H)] .= missing
        H[.!ismissing.(H) .&& isnan.(H)] .= missing
        # ensure that there are no negative widths
        W[.!ismissing.(W) .&& W .< 0] .= missing
        nid = nodes["node_id"][:]
        dmap = Dict(nid[k] => k for k=1:length(nid))
        i = [dmap[k] for k in nids]
        H[i, :], W[i, :], S[i, :], dA, Hr, Wr, Sr
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
        g = NCDatasets.group(fd, "nodes")
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
function write_output(reachid, valid, outdir, A0, n, Qa, Qu)
    outfile = joinpath(outdir, "$(reachid)_sad.nc")
    out = Dataset(outfile, "c")
    out.attrib["valid"] = valid   # FIXME Determine what is considered valid in the context of a SAD run
    defDim(out, "nt", length(Qa))
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
    close(out)
end

"""
    main()

Main driver routine.

"""
function main()
    indir = joinpath("/mnt", "data", "input")
    outdir = joinpath("/mnt", "data", "output")

    reachfile = isempty(ARGS) ? "reaches.json" : reachfile = ARGS[1]
    reachid, swotfile, sosfile, swordfile = get_reach_files(indir, reachfile)

    nids, x = river_info(reachid, swordfile)
    H, W, S, dA, Hr, Wr, Sr = read_swot_obs(swotfile, nids)
    x, H, W, S = Sad.drop_unobserved(x, H, W, S)
    A0 = missing
    n = missing
    Qa = Array{Missing}(missing, 1, size(W,1))
    Qu = Array{Missing}(missing, 1, size(W,1))
    if all(ismissing, H) || all(ismissing, W) || all(ismissing, S) || size(H, 1) <= 1
        println("$(reachid): INVALID")
        write_output(reachid, 0, outdir, A0, n, Qa, Qu)
    else
        H, S = Sad.interpolate_negative_slopes(x, H)
        Hmin = minimum(skipmissing(H[1, :]))
        Qp, np, rp, zp = Sad.priors(sosfile, Hmin, reachid)
        if ismissing(Qp)
            println("$(reachid): INVALID, missing mean discharge")
            write_output(reachid, 0, outdir, A0, n, Qa, Qu)
        else
            try
                nens = 100 # default ensemble size
                nsamples = 1000 # default sampling size
                Qa, Qu, A0, n = Sad.estimate(x, H, W, S, dA, Qp, np, rp, zp, nens, nsamples, Hr, Wr, Sr)
                println("$(reachid): VALID")
                write_output(reachid, 1, outdir, A0, n, Qa, Qu)
            catch
                println("$(reachid): INVALID")
                write_output(reachid, 0, outdir, A0, n, Qa, Qu)
            end
        end
    end
end

main()
