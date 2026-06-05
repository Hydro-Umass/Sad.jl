#!/usr/bin/env julia
#
# Science regression tests using real SWOT data
#
# These tests verify that code changes don't degrade inference performance
# on actual river reaches. They require large data files that are not
# distributed with the package, so they are skipped if data is unavailable.
#
# Run with: julia --project=. test/swot_reachtests.jl
#
# Required data files (set SWOT_DATA or use defaults):
#   - SWOT observation files (one per reach)
#   - SWORD river network files
#   - SoS priors files
#   - SVS gauge data (for metric validation)
#

using Test
using Sad
using NCDatasets
using Statistics
using Dates
using Optim
using Distributions

# --- Configuration: data paths and reach definitions ---

const DATA_ROOT = get(ENV, "SWOT_DATA", joinpath(@__DIR__, "..", "..", "data"))

const REACHES = Dict(
    12291300071 => (
        continent  = "af",
        gage_col   = 1001,   # GRDC 1291100 in SVS
        n_valid_ts_min = 30,
        expected = (
            n_range      = (0.01, 0.10),
            r_range      = (1.0, 10.0),
            z0_below_hmin = true,
            Q_median_max = 700.0,
            pearsonr_min = 0.80,
            NSE_min      = 0.50,
            KGE_min      = 0.40,
        ),
    ),
    81270600061 => (
        continent  = "na",
        gage_col   = nothing,   # USGS 4103820 not directly in SVS column map yet
        n_valid_ts_min = 30,
        expected = (
            n_range      = (0.01, 0.10),
            r_range      = (1.0, 10.0),
            z0_below_hmin = true,
            Q_median_max = 200.0,
            pearsonr_min = -Inf,   # no gage check yet
            NSE_min      = -Inf,
            KGE_min      = -Inf,
        ),
    ),
    56395000161 => (
        continent  = "oc",
        gage_col   = nothing,   # small Oceania reach, no gage
        n_valid_ts_min = 20,
        expected = (
            n_range      = (0.005, 0.15),  # wider range due to σ_obs retry (n_lo=0.005)
            r_range      = (1.0, 15.0),   # wider range due to σ_obs retry
            z0_below_hmin = false,         # z₀ may be above hmin for small rivers
            Q_median_max = 100.0,
            pearsonr_min = -Inf,           # no gage check yet
            NSE_min      = -Inf,
            KGE_min      = -Inf,
        ),
    ),
)

function _continent(reachid::Int)
    Dict(1=>"af", 2=>"eu", 3=>"as", 4=>"as", 5=>"oc",
         6=>"sa", 7=>"na", 8=>"na", 9=>"na")[floor(Int, reachid / 1e10)]
end

# --- Helper functions ---

"""Check if required data files exist for a given continent."""
function data_available(continent::String)
    sw_file = joinpath(DATA_ROOT, "sword", "$(continent)_sword_v17.nc")
    sos_file = joinpath(DATA_ROOT, "sos", "$(continent)_sword_v17_SOS_priors.nc")
    isfile(sw_file) && isfile(sos_file)
end

"""Load SWOT reach data."""
function load_reach(reachid::Int, cfg)
    continent = cfg.continent
    sw_file   = joinpath(DATA_ROOT, "sword",   "$(continent)_sword_v17.nc")
    sos_file  = joinpath(DATA_ROOT, "sos",      "$(continent)_sword_v17_SOS_priors.nc")
    swot_file = joinpath(DATA_ROOT, "swot",      "$(reachid)_SWOT.nc")
    svs_file  = joinpath(DATA_ROOT, "svs",        "SVS_v1_0.5_trans_v17.nc")

    nids, x = river_info(reachid, sw_file)
    H, W, S, dA, Hr, Wr, Sr, time_str = read_swot_obs(swot_file, nids)
    reach = Sad.preprocess(x, H, W, S)
    p = Sad.priors(sos_file, reach.hmin, reachid; S0=mean(reach.S0.(reach.x)))

    return reach, p, time_str, svs_file
end

# Reuse swot.jl I/O functions
function read_swot_obs(ncfile::String, nids::Vector{Int})
    Dataset(ncfile) do ds
        nodes = ds.group["node"]
        reaches = ds.group["reach"]
        S = permutedims(nodes["slope2"][:, :])
        H = permutedims(nodes["wse"][:, :])
        W = permutedims(nodes["width"][:, :])
        dA = convert(Vector{Sad.FloatM}, reaches["d_x_area"][:])
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

function river_info(id::Int, swordfile::String)
    Dataset(swordfile) do fd
        g = fd.group["nodes"]
        idx = findall(g["reach_id"][:] .== id)
        nid = g["node_id"][idx]
        x = g["dist_out"][idx]
        k = findall(.!ismissing.(nid))
        x = x[k]; nid = nid[k]
        x = x .- minimum(x)
        j = sortperm(x)
        convert(Vector{Int}, nid[j]), x[j]
    end
end

"""Load SVS gauge data and match to SAD time steps."""
function load_gage(svs_file::String, reachid::Int, gage_col::Int,
                   time_str::Vector{String})
    isfile(svs_file) || return nothing, nothing
    Dataset(svs_file) do ds
        Q_var = ds["Q"]
        rid_var = ds["reach_id_v17"]
        time_raw = ds["time"]

        # Verify reach matches
        matching_col = 0
        for c in 1:size(rid_var, 2)
            rids = rid_var[:, c]
            if any(.!ismissing.(rids) .&& rids .== reachid)
                # prefer the specified column if it also matches
                if c == gage_col
                    matching_col = c; break
                elseif matching_col == 0
                    matching_col = c
                end
            end
        end
        matching_col == 0 && return nothing, nothing

        Q_station = Q_var[:, matching_col]
        time_days = [ismissing(t) ? NaN : Float64(Dates.value(Date(t) - Date(2023,1,1))) for t in time_raw]
        valid_gage = findall(.!ismissing.(Q_station) .&& .!isnan.(Float64.(Q_station)) .&& Float64.(Q_station) .> 0)
        gage_Q = Float64.(Q_station[valid_gage])
        gage_dates = [Date(2023,1,1) + Day(round(Int, time_days[k])) for k in valid_gage]

        gage_dict = Dict{Date, Float64}()
        for (gi, gd) in enumerate(gage_dates)
            gage_dict[gd] = gage_Q[gi]
        end

        sad_dates = [Date(DateTime(s)) for s in time_str]
        obs_Q, mod_Q = Float64[], Float64[]
        for sad_date in sad_dates
            best_delta = 999
            best_gq = NaN
            for (gd, gq) in gage_dict
                delta = abs((sad_date - gd).value)
                if delta < best_delta
                    best_delta = delta
                    best_gq = gq
                end
            end
            if best_delta <= 1 && !isnan(best_gq)
                push!(obs_Q, best_gq)
                # placeholder for mod_Q — filled after inference
            end
        end
        return gage_dict, sad_dates
    end
end

"""Compute validation metrics between observed and modeled Q."""
function compute_metrics(obs, mod)
    n = length(obs)
    n == 0 && return (r=NaN, NSE=NaN, KGE=NaN, nBIAS=NaN, RMSE=NaN)
    mo = mean(obs); mm = mean(mod)
    so = std(obs); sm = std(mod)
    r = cor(obs, mod)
    NSE = 1 - sum((obs .- mod).^2) / sum((obs .- mo).^2)
    α = sm / so; β = mm / mo
    KGE = 1 - sqrt((r-1)^2 + (α-1)^2 + (β-1)^2)
    nBIAS = (mm - mo) / mo
    RMSE = sqrt(mean((obs .- mod).^2))
    return (r=r, NSE=NSE, KGE=KGE, nBIAS=nBIAS, RMSE=RMSE)
end

# --- Test definitions ---

@testset "SWOT Reach Regression Tests" begin

for (reachid, cfg) in sort(collect(pairs(REACHES)); by=first)
    @testset "Reach $reachid ($(cfg.continent))" begin

        # Skip if data unavailable
        swot_file = joinpath(DATA_ROOT, "swot", "$(reachid)_SWOT.nc")
        if !data_available(cfg.continent) || !isfile(swot_file)
            @info "Skipping reach $reachid — data files not found in $DATA_ROOT"
            continue
        end

        # --- Load and preprocess ---
        @testset "Data loading" begin
            reach, p, time_str, svs_file = load_reach(reachid, cfg)

            @test reach.nx > 0
            @test reach.nt > 0
            @test sum(reach.valid) >= cfg.n_valid_ts_min
            @test reach.hmin > 0
            @test all(reach.wbf.(reach.x) .> 0)
            @test all(reach.S0.(reach.x) .> 0)
        end

        # --- Inference ---
        reach, p, time_str, svs_file = load_reach(reachid, cfg)
        res = Sad.infer(p, reach; time_str=time_str)

        @testset "Inference convergence" begin
            @test res.n_post > 0
            @test res.r_post > 0
            @test isfinite(res.z0_post)
        end

        @testset "Physical parameter ranges" begin
            exp_n = cfg.expected.n_range
            @test exp_n[1] <= res.n_post <= exp_n[2]
            exp_r = cfg.expected.r_range
            @test exp_r[1] <= res.r_post <= exp_r[2]
            # z₀ should be below minimum observed WSE (hmin)
            if cfg.expected.z0_below_hmin
                @test res.z0_post < reach.hmin
            end
        end

        @testset "Discharge plausibility" begin
            valid_Q = filter(q -> q > 0, res.Q_post)
            @test length(valid_Q) >= cfg.n_valid_ts_min
            @test median(valid_Q) <= cfg.expected.Q_median_max
            @test minimum(valid_Q) > 0
        end

        # --- Gauge validation (if gauge available) ---
        if cfg.gage_col !== nothing && isfile(svs_file)
            @testset "Gage validation" begin
                gage_dict, sad_dates = load_gage(svs_file, reachid,
                                                  cfg.gage_col, time_str)

                if gage_dict !== nothing
                    # Match SAD outputs to gage
                    obs_Q, mod_Q = Float64[], Float64[]
                    for (si, sd) in enumerate(sad_dates)
                        best_delta = 999
                        best_gq = NaN
                        for (gd, gq) in gage_dict
                            delta = abs((sd - gd).value)
                            if delta < best_delta
                                best_delta = delta
                                best_gq = gq
                            end
                        end
                        if best_delta <= 1 && !isnan(best_gq) && res.Q_post[si] > 0
                            push!(obs_Q, best_gq)
                            push!(mod_Q, res.Q_post[si])
                        end
                    end

                    if length(obs_Q) >= 10
                        metrics = compute_metrics(obs_Q, mod_Q)
                        @info "Reach $reachid metrics" metrics...
                        @test metrics.r >= cfg.expected.pearsonr_min
                        @test metrics.NSE >= cfg.expected.NSE_min
                        @test metrics.KGE >= cfg.expected.KGE_min
                    else
                        @warn "Reach $reachid: only $(length(obs_Q)) matched gage timesteps (need ≥10)"
                    end
                else
                    @warn "Reach $reachid: no gage data found at column $(cfg.gage_col)"
                end
            end
        end

        # --- σ_obs estimation ---
        @testset "σ_obs estimation" begin
            months = map(time_str) do s
                try Month(DateTime(s)).value; catch; 1; end
            end
            precomp = Sad.precompute_manning(reach, months;
                S0=mean(reach.S0.(reach.x)))
            @test precomp.σ_obs_est > 0
            @test precomp.σ_obs_est >= 0.3  # floor value
        end

        # --- Prior sanity ---
        @testset "Prior sanity" begin
            # Q priors should have positive support
            for mo in 1:12
                @test quantile(p.Qp[mo], 0.01) > 0
            end
            # n prior should be in physically plausible range
            @test minimum(p.np) >= 0.01
            @test maximum(p.np) <= 0.2
            # r prior should allow positive values
            @test minimum(p.rp) > 0
            # z₀ prior should have support below hmin
            @test maximum(p.zp) <= reach.hmin
        end
    end
end

end # testset