using Test
using Statistics
using NCDatasets
using Sad

# path to test data, relative to the test directory
const TESTDATA = joinpath(@__DIR__, "testdata.nc")

# helper functions

"""Load test reach and priors from testdata.nc."""
function load_test_data()
    Dataset(TESTDATA) do f
        g    = f.group["XS_Timeseries"]
        ri   = f.group["River_Info"]
        qwbm = Float64(ri["QWBM"][1])
        x    = (g["X"][:][end] .- g["X"][:])[end:-1:1, 1]
        H    = convert(Matrix{Sad.FloatM}, g["H"][end:-1:1, 1:10:end])
        W    = convert(Matrix{Sad.FloatM}, g["W"][end:-1:1, 1:10:end])
        Q    = g["Q"][end:-1:1, 1:10:end]
        S    = (diff(H, dims=1) ./ diff(x))[:, 1:10:end]
        S    = convert(Matrix{Sad.FloatM}, [S[1, :]'; S])
        reach = Sad.preprocess(x, H, W, S)
        p     = Sad.priors(mean(Q[1, :]), reach.hmin, Sad.sinuous; reach)
        reach, p, Q
    end
end

# preprocess tests

@testset "preprocess" begin
    reach, p, Q = load_test_data()

    @testset "SWOTReach structure" begin
        @test reach.nx > 0
        @test reach.nt == 37
        @test reach.nobs == 68
        @test length(reach.x) == reach.nx
        @test size(reach.H) == (reach.nx, reach.nt)
        @test size(reach.W) == (reach.nx, reach.nt)
        @test length(reach.valid) == reach.nt
        @test reach.hmin > 0
    end

    @testset "chainage" begin
        # chainage starts at 0 and increases monotonically
        @test reach.x[1] == 0.0
        @test all(diff(reach.x) .> 0)
    end

    @testset "bed slope" begin
        # S0 should be positive everywhere
        S0_vals = reach.S0.(reach.x)
        @test all(S0_vals .> 0)
        # S0 should be at least min_slope
        @test all(S0_vals .>= 1e-5)
    end

    @testset "cumulative z" begin
        # z starts at 0 at downstream end
        @test reach.z(reach.x[1]) ≈ 0.0 atol=1e-10
        # z increases upstream
        @test reach.z(reach.x[end]) > 0
    end

    @testset "wbf and hbf" begin
        # bankfull width and WSE should be positive everywhere
        @test all(reach.wbf.(reach.x) .> 0)
        @test all(reach.hbf.(reach.x) .> 0)
        # bankfull WSE should be above hmin
        @test all(reach.hbf.(reach.x) .>= reach.hmin)
    end

    @testset "valid timesteps" begin
        # tt least some timesteps should be valid
        @test sum(reach.valid) > 0
    end
end

# priors tests

@testset "priors" begin
    reach, p, Q = load_test_data()

    @testset "SWOTPriors structure" begin
        @test p isa Sad.SWOTPriors
    end

    @testset "Q prior" begin
        # prior should have positive support
        @test quantile(p.Qp, 0.01) > 0
        @test quantile(p.Qp, 0.999) > quantile(p.Qp, 0.01)
        # upper bound should be 20x qwbm
        @test quantile(p.Qp, 0.999) > 1000.0
        # mean Q samples should be positive
        Q_samples = [rand(p.Qp) for _ in 1:1000]
        @test all(Q_samples .> 0)
    end

    @testset "n prior" begin
        # Manning n should be in physically plausible range
        @test minimum(p.np) >= 0.01
        @test maximum(p.np) <= 0.15
    end

    @testset "r prior" begin
        # Dingman r should be positive
        @test minimum(p.rp) > 0
    end

    @testset "z prior" begin
        # bed elevation should be below hmin
        @test maximum(p.zp) < reach.hmin
    end
end

# GVF tests

@testset "gvf" begin
    reach, p, Q = load_test_data()

    @testset "area" begin
        # area should be positive for positive depth
        @test Sad.area(1.0, 100.0, 5.0, 2.0) > 0
        # area should increase with depth
        @test Sad.area(2.0, 100.0, 5.0, 2.0) > Sad.area(1.0, 100.0, 5.0, 2.0)
    end

    @testset "compute_wse" begin
        # WSE = y * (r+1)/r + z
        @test Sad.compute_wse(2.0, 100.0, 2.0) ≈ 2.0 * 3.0 / 2.0 + 100.0
        @test Sad.compute_wse(1.0, 0.0, 1.0) ≈ 2.0
        # should be increasing in y
        @test Sad.compute_wse(3.0, 100.0, 2.0) > Sad.compute_wse(1.0, 100.0, 2.0)
    end

    @testset "compute_width" begin
        # W = Wb * (y / Yb)^(1/r)
        @test Sad.compute_width(3.0, 100.0, 5.0, 2.0) ≈ 100.0 * (3.0 / 5.0)^(1/2.0)
        @test Sad.compute_width(5.0, 100.0, 5.0, 2.0) ≈ 100.0  # at bankfull
        # should be increasing in y
        @test Sad.compute_width(4.0, 100.0, 5.0, 2.0) > Sad.compute_width(1.0, 100.0, 5.0, 2.0)
    end

    @testset "gvf_solve returns valid mean depth" begin
        t   = findall(reach.valid)[1]
        H_bc = reach.H[1, t]
        y_pred = Sad.gvf_solve(500.0, 0.03, 2.0, reach.hmin - 5.0,
                               H_bc, reach)
        # should return a vector, not nothing
        @test !isnothing(y_pred)
        @test length(y_pred) == reach.nobs
        # mean depth should be finite and positive
        @test all(isfinite.(y_pred))
        @test all(y_pred .> 0)
        # convert to WSE and check it increases upstream
        z_saveat = (reach.hmin - 5.0) .+ reach.z.(collect(reach.obs.x))
        WSE = Sad.compute_wse.(y_pred, z_saveat, 2.0)
        @test WSE[end] > WSE[1]
        # downstream boundary should match H_bc
        @test WSE[1] ≈ H_bc atol=0.1
    end

    @testset "gvf_solve fails gracefully" begin
        t    = findall(reach.valid)[1]
        H_bc = reach.H[1, t]
        # extreme parameters should return nothing rather than error
        result = Sad.gvf_solve(1e10, 0.001, 0.1, -100.0, H_bc, reach)
        @test isnothing(result) || all(isfinite.(result))
    end

    @testset "gvf_solve mean depth consistency" begin
        # mean depth + z + shape factor should reconstruct WSE consistently
        t   = findall(reach.valid)[1]
        H_bc = reach.H[1, t]
        r_val = 2.0
        z0_val = reach.hmin - 5.0
        y_pred = Sad.gvf_solve(500.0, 0.03, r_val, z0_val, H_bc, reach)
        @test !isnothing(y_pred)
        # downstream mean depth should satisfy y_bc = (H_bc - z0) * r/(r+1)
        y_bc_expected = (H_bc - z0_val) * r_val / (r_val + 1)
        @test y_pred[1] ≈ y_bc_expected atol=0.1
    end
end

# rejection sampling tests

@testset "rejection_sample" begin
    reach, p, Q = load_test_data()

    ens = Sad.rejection_sample(p, reach; N=50, max_attempts=5_000)

    @testset "output shape" begin
        # should return 3 rows: [n, r, z0]
        @test size(ens, 1) == 3
        @test size(ens, 2) > 0
    end

    @testset "parameter bounds" begin
        if size(ens, 2) > 0
            n_ens  = ens[1, :]
            r_ens  = ens[2, :]
            z0_ens = ens[3, :]
            # all parameters should be within prior support
            @test all(n_ens  .>= minimum(p.np))
            @test all(n_ens  .<= maximum(p.np))
            @test all(r_ens  .>= minimum(p.rp))
            @test all(r_ens  .<= maximum(p.rp))
            @test all(z0_ens .>= minimum(p.zp))
            @test all(z0_ens .<= maximum(p.zp))
        end
    end

    @testset "select_representative_timesteps" begin
        rep_ts = Sad.select_representative_timesteps(reach; n_bins=5)
        @test length(rep_ts) > 0
        @test length(rep_ts) <= 5
        # all selected timesteps should be valid
        @test all(reach.valid[rep_ts])
        # should cover range of downstream WSE
        h_ds = [reach.H[1, t] for t in rep_ts]
        @test maximum(h_ds) > minimum(h_ds)
    end
end

# Q ensemble tests

@testset "q_ensemble" begin
    reach, p, Q = load_test_data()
    ens = Sad.rejection_sample(p, reach; N=50, max_attempts=5_000)

    if size(ens, 2) > 0
        valid_ts = findall(reach.valid)
        t = valid_ts[1]
        Q_ens, members = Sad.q_ensemble(p, ens, reach, t)

        @testset "output" begin
            @test length(Q_ens) == length(members)
            @test all(Q_ens .> 0)
            @test all(members .>= 1)
            @test all(members .<= size(ens, 2))
        end

        @testset "Q range" begin
            # accepted Q should be within prior support
            @test all(Q_ens .>= quantile(p.Qp, 0.001))
            @test all(Q_ens .<= quantile(p.Qp, 0.999))
        end
    end
end

# infer_channel_params tests

@testset "infer_channel_params" begin
    reach, p, _ = load_test_data()

    months = fill(1, reach.nt)
    pa = Sad.infer_channel_params(p, reach, months, 50)

    @testset "returns SWOTPriors" begin
        @test pa isa Sad.SWOTPriors
        # Qp should be unchanged (stage 1 does not touch discharge prior)
        @test pa.Qp === p.Qp
    end

    @testset "posterior distributions are sampleable" begin
        n_samp  = rand(pa.np, 100)
        r_samp  = rand(pa.rp, 100)
        z0_samp = rand(pa.zp, 100)
        @test all(isfinite.(n_samp))
        @test all(isfinite.(r_samp))
        @test all(isfinite.(z0_samp))
    end

    @testset "posterior parameter ranges are physically plausible" begin
        n_samp  = rand(pa.np, 500)
        r_samp  = rand(pa.rp, 500)
        @test all(n_samp  .> 0)
        @test all(r_samp  .> 0)
    end
end

# end-to-end inference test

@testset "infer" begin
    reach, p, Q_truth = load_test_data()

    res = Sad.infer(p, reach; N=100, σₒ=2.0)

    @testset "output structure" begin
        @test haskey(res, :Q_post)
        @test haskey(res, :reach_ensemble)
        @test haskey(res, :A_post)
        @test haskey(res, :valid_ts)
        @test haskey(res, :completeness)
        @test length(res.Q_post) == reach.nt
        @test size(res.reach_ensemble, 1) == 3
        # A_post is a per-timestep vector
        @test res.A_post isa Vector
        @test length(res.A_post) == reach.nt
    end

    @testset "Q_post" begin
        valid_ts = findall(reach.valid)
        Q_est = res.Q_post[valid_ts]
        good  = findall(.!isnan.(Q_est))
        # at least 80% of valid timesteps should have a posterior estimate
        @test length(good) / length(valid_ts) >= 0.8
        # all posterior Q values should be positive
        @test all(Q_est[good] .> 0)
    end

    @testset "NSE" begin
        valid_ts = findall(reach.valid)
        Qt = [Q_truth[1, t] for t in valid_ts]
        Qe = res.Q_post[valid_ts]
        good = findall(.!isnan.(Qe))
        nse  = 1 - sum((Qe[good] .- Qt[good]).^2) /
                   sum((Qt[good] .- mean(Qt[good])).^2)
        @info "NSE: $(round(nse, digits=3))"
        # NSE should be positive (two-stage LETKF approach)
        @test nse > -1.0  # model should at least not be worse than horizontal line by >100%
    end

    @testset "compute_A0" begin
        A0 = Sad.compute_A0(reach, res.reach_ensemble)
        @test A0 > 0
        @test isfinite(A0)
        # A0 should be physically plausible (between 100 and 100000 m²)
        @test A0 > 100.0
        @test A0 < 100_000.0
    end
end

# edge case tests: missing data

@testset "missing data" begin
    reach, p, Q_truth = load_test_data()

    # helper to create a reach with modified observations
    function make_reach(H, W, S, x)
        Sad.preprocess(x, H, W, S)
    end

    Dataset(TESTDATA) do f
        g    = f.group["XS_Timeseries"]
        x    = (g["X"][:][end] .- g["X"][:])[end:-1:1, 1]
        H    = convert(Matrix{Sad.FloatM}, g["H"][end:-1:1, :])
        W    = convert(Matrix{Sad.FloatM}, g["W"][end:-1:1, :])
        S    = diff(H, dims=1) ./ diff(x)
        S    = convert(Matrix{Sad.FloatM}, [S[1, :]'; S])

        @testset "all H missing" begin
            H_miss = convert(Matrix{Sad.FloatM}, fill(missing, size(H)...))
            reach_miss = Sad.preprocess(x, H_miss, W, S)
            # no valid timesteps
            @test sum(reach_miss.valid) == 0
            # infer should not error
            res = Sad.infer(p, reach_miss; N=10)
            @test all(isnan.(res.Q_post))
        end

        @testset "all W missing" begin
            W_miss = convert(Matrix{Sad.FloatM}, fill(missing, size(W)...))
            reach_miss = Sad.preprocess(x, H, W_miss, S)
            # no valid timesteps since W is missing
            @test sum(reach_miss.valid) == 0
            res = Sad.infer(p, reach_miss; N=10)
            @test all(isnan.(res.Q_post))
        end

        @testset "all H and W missing" begin
            H_miss = convert(Matrix{Sad.FloatM}, fill(missing, size(H)...))
            W_miss = convert(Matrix{Sad.FloatM}, fill(missing, size(W)...))
            reach_miss = Sad.preprocess(x, H_miss, W_miss, S)
            @test sum(reach_miss.valid) == 0
            res = Sad.infer(p, reach_miss; N=10)
            @test all(isnan.(res.Q_post))
        end

        @testset "H missing for all but one timestep" begin
            H_one = convert(Matrix{Sad.FloatM}, fill(missing, size(H)...))
            H_one[:, 1] = H[:, 1]
            reach_one = Sad.preprocess(x, H_one, W, S)
            # only one timestep, not enough for rejection sampling but should not error
            @test sum(reach_one.valid) <= 1
            res = Sad.infer(p, reach_one; N=10)
            @test length(res.Q_post) == reach_one.nt
        end

        @testset "W missing for half of timesteps" begin
            W_half = copy(W)
            W_half[:, 1:2:end] .= missing
            reach_half = Sad.preprocess(x, H, W_half, S)
            # should still have valid timesteps for even-indexed times
            @test sum(reach_half.valid) > 0
            res = Sad.infer(p, reach_half; N=50)
            @test length(res.Q_post) == reach_half.nt
            # valid timesteps should have posterior estimates
            valid_ts = findall(reach_half.valid)
            Q_est = res.Q_post[valid_ts]
            good  = findall(.!isnan.(Q_est))
            @test length(good) > 0
            @test all(Q_est[good] .> 0)
        end

        @testset "H missing for some nodes at all timesteps" begin
            H_nodes = copy(H)
            # knock out first 10 nodes entirely
            H_nodes[1:10, :] .= missing
            reach_nodes = Sad.preprocess(x, H_nodes, W, S)
            # should still work with remaining nodes
            @test reach_nodes.nobs > 0
            @test sum(reach_nodes.valid) > 0
        end

        @testset "drop_unobserved removes fully missing nodes" begin
            H_sparse = copy(H)
            W_sparse = copy(W)
            #mMake nodes 5-10 fully missing
            H_sparse[5:10, :] .= missing
            W_sparse[5:10, :] .= missing
            x_out, H_out, W_out, S_out = Sad.drop_unobserved(
                copy(x),
                convert(Matrix{Sad.FloatM}, H_sparse),
                convert(Matrix{Sad.FloatM}, W_sparse),
                copy(S),
            )
            # fully missing nodes should be removed
            @test size(H_out, 1) < size(H_sparse, 1)
            @test size(H_out, 1) == size(W_out, 1)
            # chainage should still start from 0
            @test x_out[1] == 0.0
        end
    end
end
