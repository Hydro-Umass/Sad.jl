### Run SAD algorithm with SWOT data

using Sad
using ArgParse
using DelimitedFiles
using Distributions
using JSON
using LinearAlgebra
using NCDatasets
using PyCall

const FILL = -999999999999

"""
    get_reach_files(indir, reachjson)

Get reach file names and download SoS file.

"""
function get_reach_files(indir, tmpdir, reachjson, index, sosbucket)
    data = Ref{Dict{String, Any}}
    open(joinpath(indir, reachjson)) do jf
        data = read(jf, String)
    end

    reachlist = JSON.parse(data)[index]
    if !isempty(sosbucket)
        sosfile = joinpath(tmpdir, reachlist["sos"])
        pushfirst!(pyimport("sys")."path", "/app/sos_read")   # Load sos_read script
        downloadsos = pyimport("sos_read")["download_sos"]
        downloadsos(sosbucket, sosfile)
        reachlist["reach_id"], joinpath(indir, "swot", reachlist["swot"]), sosfile, joinpath(indir, "sword", reachlist["sword"])
    else
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
        nid = nodes["node_id"][:]
        dmap = Dict(nid[k] => k for k=1:length(nid))
        i = [dmap[k] for k in nids]
        time_str_var = reaches["time_str"].var
        time_str_raw = permutedims(time_str_var[:])
        time_str = [join(time_str_raw[i, :]) for i in 1:size(time_str_raw, 1)]


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
    parse_commandline(r)

Parse command line for arguments.

"""
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--index", "-i"
            help = "Index of reach to run on"
            arg_type = Int
            default = 0
        "--reachfile", "-r"
            help = "Name of reaches JSON file"
            arg_type = String
            default = "reaches.json"
        "--bucketkey", "-b"
            help = "Bucket and key prefix to download SoS from"
            arg_type = String
            default = ""
    end

    return parse_args(s)
end

"""
    main()

Main driver routine.

"""
function main()
    indir = joinpath("/mnt", "data", "input")
    outdir = joinpath("/mnt", "data", "output")
    tmpdir = joinpath("/tmp")

    parsed_args = parse_commandline()
    if parsed_args["index"] == -256
        index = try parse(Int64, ENV["AWS_BATCH_JOB_ARRAY_INDEX"]) + 1 catch KeyError 1 end
    else
        index = parsed_args["index"] + 1
    end

    reachfile = parsed_args["reachfile"]
    bucketkey = parsed_args["bucketkey"]
    println("Index: $(index)")
    println("Reach File: $(reachfile)")
    println("Bucket Key: $(bucketkey)")

    reachid, swotfile, sosfile, swordfile = get_reach_files(indir, tmpdir, reachfile, index, bucketkey)
    println("Reach ID: $(reachid)")
    println("SWOT: $(swotfile)")
    println("SOS: $(sosfile)")
    println("SWORD: $(swordfile)")

    nids, x = river_info(reachid, swordfile)
    H, W, S, dA, Hr, Wr, Sr, time_str = read_swot_obs(swotfile, nids)

    try
        reach = Sad.preprocess(x, H, W, S)
    catch e
        if e isa MethodError
            println("Error loading swot observation")
            end
        end
    A0 = missing
    n = missing
    Qa = Array{Missing}(missing, 1, size(W, 2))
    Qu = Array{Missing}(missing, 1, size(W, 2))
    if all(ismissing, H) || all(ismissing, W) || all(ismissing, S)
        println("$(reachid): INVALID")
        write_output(reachid, 0, outdir, A0, n, Qa, Qu, W, time_str)
    else
        p = Sad.priors(sosfile, reach.hmin, reachid)
        if ismissing(p.Qp)
            println("$(reachid): INVALID, missing mean discharge")
            write_output(reachid, 0, outdir, A0, n, Qa, Qu, W, time_str)
        else
            try
                res = Sad.infer(p, reach)
                A0  = Sad.compute_A0(reach, res.reach_ensemble)
                n   = mean(res.reach_ensemble[1, :])
                Qa[1, :]  = res.Q_post
                Qu[1, :]  = [isnothing(res.A_post[t]) ? NaN : std(exp.(res.A_post[t][1,:])) for t in 1:reach.nt]
                println("$(reachid): VALID")
                write_output(reachid, 1, outdir, A0, n, Qa, Qu, W, time_str)
            catch
                println("$(reachid): INVALID")
                write_output(reachid, 0, outdir, A0, n, Qa, Qu, W, time_str)
            end
        end
    end
end

main()
