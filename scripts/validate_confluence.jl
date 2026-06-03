using NCDatasets, Dates, DataFrames

include("swot.jl")

kge(o, m) = 1 - sqrt((cor(o, m) - 1)^2 + (std(m) / std(o) - 1)^2 + (mean(m) / mean(o) - 1)^2)

nse(o, m) = 1 - sum((o .- m).^2) / sum((o .- mean(o)).^2)

continent(reachid::Integer) = Dict(
    1 => "af", 2 => "eu", 3 => "as", 4 => "as", 5 => "oc",
    6 => "sa", 7 => "na", 8 => "na", 9 => "na",
)[floor(Int, reachid / 1e10)]

function test_confluence(reachid)
    swordfile = "../../data/sword/$(continent(reachid))_sword_v17.nc"
    sosfile = "../../data/sos/$(continent(reachid))_sword_v17_SOS_priors.nc"
    swotfile = "../../data/swot/$(reachid)_SWOT.nc"
    main(reachid, swordfile, sosfile, swotfile, ".")
    Qa, ta = Dataset("./$(reachid)_sad.nc") do ds
        ds["Qa"][:], [t == "missing" ? missing : Date(DateTime(t)) for t in ds["time_str"][:]]
    end
    svsfile = "../../data/svs/SVS_v1_0.5_trans_v17.nc"
    Qobs, tobs = Dataset(svsfile) do ds
        rids = ds["reach_id_v17"][1, :]
        i = findfirst(x -> !ismissing(x) && x==reachid, rids)
        ds["Q"][:, i], Date.(ds["time"][:])
    end
    return DataFrame(Q=Qa, t=ta) |> dropmissing, DataFrame(Q=Qobs, t=tobs) |> dropmissing
end

function get_v1(reachid)
    Q1, t1 = Dataset("../../outputv1/$(reachid)_sad.nc") do ds
        ds["Qa"][:], [t == "missing" ? missing : Date(DateTime(t)) for t in ds["time_str"][:]]
    end
    return DataFrame(Q=Q1, t=t1) |> dropmissing
end
