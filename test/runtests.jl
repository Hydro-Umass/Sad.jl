using Test
using Statistics
using LinearAlgebra
using NCDatasets
using ForwardDiff
using Distributions
using BSplineKit
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
        S_row = S[1:1, :]
        S    = convert(Matrix{Sad.FloatM}, vcat(S_row, S))
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
        # at least some timesteps should be valid
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

# Dingman cross-section geometry tests

@testset "cross-section geometry" begin
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
end

# Manning analytical model tests (v2)

@testset "manning" begin
    reach, p, Q_truth = load_test_data()

    @testset "manning_Q" begin
        # basic sanity: positive Q for positive inputs
        @test Sad.manning_Q(0.03, 2.0, 5.0, 200.0, 10.0, 1e-4) > 0
        # Q increases with depth
        @test Sad.manning_Q(0.03, 2.0, 6.0, 200.0, 10.0, 1e-4) >
              Sad.manning_Q(0.03, 2.0, 5.0, 200.0, 10.0, 1e-4)
        # Q decreases with roughness
        @test Sad.manning_Q(0.04, 2.0, 5.0, 200.0, 10.0, 1e-4) <
              Sad.manning_Q(0.03, 2.0, 5.0, 200.0, 10.0, 1e-4)
        # Q increases with slope
        @test Sad.manning_Q(0.03, 2.0, 5.0, 200.0, 10.0, 2e-4) >
              Sad.manning_Q(0.03, 2.0, 5.0, 200.0, 10.0, 1e-4)
        # Q increases with width
        @test Sad.manning_Q(0.03, 2.0, 5.0, 300.0, 10.0, 1e-4) >
              Sad.manning_Q(0.03, 2.0, 5.0, 200.0, 10.0, 1e-4)
        # consistency with Dingman area function
        n, r, y, Wb, Yb, S = 0.03, 2.0, 5.0, 200.0, 10.0, 1e-4
        A = Sad.area(y, Wb, Yb, r)
        Q_manning = (1/n) * A * y^(2/3) * sqrt(S)
        @test Sad.manning_Q(n, r, y, Wb, Yb, S) ≈ Q_manning
    end

    @testset "manning_depth" begin
        # inverses manning_Q
        n, r, Wb, Yb, S = 0.03, 2.0, 200.0, 10.0, 1e-4
        for Q in [100.0, 500.0, 2000.0]
            y = Sad.manning_depth(Q, n, r, Wb, Yb, S)
            @test y > 0
            @test Sad.manning_Q(n, r, y, Wb, Yb, S) ≈ Q
        end
        # depth increases with Q
        @test Sad.manning_depth(1000.0, 0.03, 2.0, 200.0, 10.0, 1e-4) >
              Sad.manning_depth(500.0, 0.03, 2.0, 200.0, 10.0, 1e-4)
        # depth decreases with slope
        @test Sad.manning_depth(500.0, 0.03, 2.0, 200.0, 10.0, 2e-4) <
              Sad.manning_depth(500.0, 0.03, 2.0, 200.0, 10.0, 1e-4)
    end

    @testset "backwater_length" begin
        # positive for positive inputs
        @test Sad.backwater_length(5.0, 1e-4, 2.0) > 0
        # increases with depth (deeper → longer backwater)
        @test Sad.backwater_length(10.0, 1e-4, 2.0) >
              Sad.backwater_length(5.0, 1e-4, 2.0)
        # decreases with slope (steeper → shorter backwater)
        @test Sad.backwater_length(5.0, 2e-4, 2.0) <
              Sad.backwater_length(5.0, 1e-4, 2.0)
        # rectangular limit: λ → 3*y0/(10*S0) as r → ∞
        y0, S0 = 5.0, 1e-4
        λ_rect = 3 * y0 / (10 * S0)
        λ_r100 = Sad.backwater_length(y0, S0, 100.0)
        @test λ_r100 ≈ λ_rect rtol=0.01  # within 1% for r=100
        # Dingman correction: λ_r < λ_rect for finite r
        @test Sad.backwater_length(y0, S0, 2.0) < λ_rect
    end

    @testset "manning_wse" begin
        # uniform flow WSE should increase upstream (bed rises)
        n, r, z0, S = 0.03, 2.0, 100.0, 1e-4
        x_nodes = [0.0, 2000.0, 4000.0, 6000.0]
        Wb_nodes = fill(200.0, 4)
        Yb_nodes = fill(10.0, 4)
        z_nodes = [0.0, 0.2, 0.4, 0.6]
        wse = Sad.manning_wse(500.0, n, r, z0, S, x_nodes, Wb_nodes, Yb_nodes, z_nodes)
        @test length(wse) == 4
        @test all(isfinite.(wse))
        # WSE increases upstream due to bed rise
        @test all(diff(wse) .> 0)
        # at x=0: WSE = z0 + 0 + y0*(r+1)/r
        y0 = Sad.manning_depth(500.0, n, r, 200.0, 10.0, S)
        @test wse[1] ≈ z0 + y0 * (r + 1) / r
    end

    @testset "manning_wse_backwater boundary conditions" begin
        n, r, z0, S0 = 0.03, 2.0, 100.0, 1e-4
        x_nodes = [0.0, 2000.0, 4000.0, 6000.0]
        Wb_nodes = fill(200.0, 4)
        Yb_nodes = fill(10.0, 4)
        z_nodes = [0.0, 0.2, 0.4, 0.6]
        y0 = Sad.manning_depth(500.0, n, r, 200.0, 10.0, S0)

        # H_bc = uniform flow WSE → backwater = uniform flow
        H_bc_unif = Sad.compute_wse(y0, z0, r)
        wse_unif = Sad.manning_wse(500.0, n, r, z0, S0, x_nodes, Wb_nodes, Yb_nodes, z_nodes)
        wse_bw = Sad.manning_wse_backwater(500.0, n, r, z0, S0, H_bc_unif,
                                            x_nodes, Wb_nodes, Yb_nodes, z_nodes)
        @test wse_bw ≈ wse_unif atol=1e-10

        # H_bc above uniform → backwater: WSE[1] = H_bc exactly
        H_bc_high = H_bc_unif + 2.0
        wse_bw2 = Sad.manning_wse_backwater(500.0, n, r, z0, S0, H_bc_high,
                                             x_nodes, Wb_nodes, Yb_nodes, z_nodes)
        @test wse_bw2[1] ≈ H_bc_high
        # upstream WSE should approach uniform flow
        @test wse_bw2[end] > wse_unif[end]
        @test wse_bw2[end] < H_bc_high + z_nodes[end]

        # H_bc below uniform → drawdown (M2): WSE[1] = H_bc exactly
        H_bc_low = H_bc_unif - 1.0
        wse_bw3 = Sad.manning_wse_backwater(500.0, n, r, z0, S0, H_bc_low,
                                             x_nodes, Wb_nodes, Yb_nodes, z_nodes)
        @test wse_bw3[1] ≈ H_bc_low
    end

    @testset "ForwardDiff compatibility" begin
        # manning_Q and manning_depth must work with Dual numbers
        n, r, Wb, Yb, S = 0.03, 2.0, 200.0, 10.0, 1e-4
        Q = 500.0

        # scalar derivatives
        dQ_dy = ForwardDiff.derivative(y -> Sad.manning_Q(n, r, y, Wb, Yb, S), 5.0)
        @test isfinite(dQ_dy)
        @test dQ_dy > 0  # Q increases with y

        dy_dQ = ForwardDiff.derivative(Q -> Sad.manning_depth(Q, n, r, Wb, Yb, S), Q)
        @test isfinite(dy_dQ)
        @test dy_dQ > 0  # y increases with Q

        # manning_wse_backwater gradient w.r.t. Q
        x_nodes = [0.0, 2000.0, 4000.0]
        Wb_nodes = [200.0, 200.0, 200.0]
        Yb_nodes = [10.0, 10.0, 10.0]
        z_nodes = [0.0, 0.2, 0.4]
        z0, S0, H_bc = 100.0, 1e-4, 106.0

        grad = ForwardDiff.derivative(
            Q -> Sad.manning_wse_backwater(Q, n, r, z0, S0, H_bc,
                                           x_nodes, Wb_nodes, Yb_nodes, z_nodes)[2],
            Q)
        @test isfinite(grad)

        # multi-parameter gradient (as Optim.L-BFGS would use)
        function obj(theta)
            wse = Sad.manning_wse_backwater(
                theta[1], n, r, theta[2], S0, H_bc,
                x_nodes, Wb_nodes, Yb_nodes, z_nodes)
            return sum((wse[2:end] .- 107.0).^2)
        end
        g = ForwardDiff.gradient(obj, [Q, z0])
        @test all(isfinite.(g))
        @test length(g) == 2
    end
end

# v2 inference tests

@testset "infer" begin
    reach, p, Q_truth = load_test_data()
    valid_ts = findall(reach.valid)
    Qt = [Q_truth[1, t] for t in valid_ts]

    @testset "precompute_manning" begin
        months = fill(1, reach.nt)
        precomp = Sad.precompute_manning(reach, months)
        @test length(precomp.x_nodes) == reach.nobs
        @test precomp.nt == reach.nt
        @test length(precomp.valid_ts) == sum(reach.valid)
        @test length(precomp.H_bc) == sum(reach.valid)
        @test all(precomp.H_bc .> 0)
        @test precomp.S0 > 0
        # upstream indices should be >= 2 (j=1 is boundary)
        for js in precomp.upstream_j
            @test all(js .>= 2)
        end
    end

    @testset "neg_log_joint" begin
        months = fill(1, reach.nt)
        precomp = Sad.precompute_manning(reach, months)
        nt = reach.nt
        theta0 = zeros(nt + 3)
        for t in 1:nt
            theta0[t] = log(quantile(p.Qp, 0.5))
        end
        theta0[nt+1] = log(0.03)
        theta0[nt+2] = log(2.0)
        theta0[nt+3] = mean(p.zp)

        nll = Sad.neg_log_joint(theta0, precomp, p)
        @test isfinite(nll)
        @test nll > 0

        g = ForwardDiff.gradient(theta -> Sad.neg_log_joint(theta, precomp, p), theta0)
        @test length(g) == nt + 3
        @test all(isfinite.(g))
    end

    @testset "infer basic" begin
        res = Sad.infer(p, reach; σ_obs=1.0, ν=3.0, λ_smooth=1.0,
                            iterations=100, g_tol=1e-6)

        @test length(res.Q_post) == reach.nt
        @test length(res.Q_std) == reach.nt
        @test length(res.θ_map) == reach.nt + 3
        @test size(res.Σ) == (reach.nt + 3, reach.nt + 3)

        Qe = res.Q_post[valid_ts]
        @test all(Qe .> 0)

        @test res.n_post > 0
        @test res.r_post > 0
        @test isfinite(res.z0_post)
    end

    @testset "infer NSE" begin
        res = Sad.infer(p, reach; σ_obs=1.0, ν=3.0, λ_smooth=5.0,
                            iterations=200)
        Qe = res.Q_post[valid_ts]
        nse = 1 - sum((Qe .- Qt).^2) / sum((Qt .- mean(Qt)).^2)
        @info "v2 NSE" nse=round(nse, digits=3)
        @test nse > -2.0  # model should be reasonable
    end

    @testset "infer speed" begin
        t_v2 = @elapsed Sad.infer(p, reach; σ_obs=1.0, ν=3.0, λ_smooth=1.0,
                                       iterations=100)
        @info "v2 inference time" t_v2=round(t_v2, digits=2)
        @test t_v2 < 30.0  # should complete in under 30 seconds
    end

    @testset "safe_logpdf" begin
        d = Uniform(5.0, 11.0)
        @test isfinite(Sad.safe_logpdf(d, 4.0))
        @test isfinite(Sad.safe_logpdf(d, 12.0))
        @test isfinite(Sad.safe_logpdf(d, 8.0))

        d_tn = truncated(Normal(0.03, 0.01), 0.01, 0.05)
        @test isfinite(Sad.safe_logpdf(d_tn, 0.06))
    end
end

# edge case tests: missing data

@testset "missing data" begin
    reach, p, Q_truth = load_test_data()

    Dataset(TESTDATA) do f
        g    = f.group["XS_Timeseries"]
        x    = (g["X"][:][end] .- g["X"][:])[end:-1:1, 1]
        H    = convert(Matrix{Sad.FloatM}, g["H"][end:-1:1, :])
        W    = convert(Matrix{Sad.FloatM}, g["W"][end:-1:1, :])
        S    = (diff(H, dims=1) ./ diff(x))
        S_row = S[1:1, :]
        S    = convert(Matrix{Sad.FloatM}, vcat(S_row, S))

        @testset "all H missing" begin
            H_miss = convert(Matrix{Sad.FloatM}, fill(missing, size(H)...))
            reach_miss = Sad.preprocess(x, H_miss, W, S)
            # no valid timesteps
            @test sum(reach_miss.valid) == 0
            # infer should not error
            res = Sad.infer(p, reach_miss; iterations=10)
            @test all(isnan.(res.Q_post))
        end

        @testset "all W missing" begin
            W_miss = convert(Matrix{Sad.FloatM}, fill(missing, size(W)...))
            reach_miss = Sad.preprocess(x, H, W_miss, S)
            # no valid timesteps since W is missing
            @test sum(reach_miss.valid) == 0
            res = Sad.infer(p, reach_miss; iterations=10)
            @test all(isnan.(res.Q_post))
        end

        @testset "all H and W missing" begin
            H_miss = convert(Matrix{Sad.FloatM}, fill(missing, size(H)...))
            W_miss = convert(Matrix{Sad.FloatM}, fill(missing, size(W)...))
            reach_miss = Sad.preprocess(x, H_miss, W_miss, S)
            @test sum(reach_miss.valid) == 0
            res = Sad.infer(p, reach_miss; iterations=10)
            @test all(isnan.(res.Q_post))
        end

        @testset "drop_unobserved removes fully missing nodes" begin
            H_sparse = copy(H)
            W_sparse = copy(W)
            # make nodes 5-10 fully missing
            H_sparse[5:10, :] .= missing
            W_sparse[5:10, :] .= missing
            x_out, H_out, W_out, S_out = Sad.drop_unobserved(
                copy(x),
                convert(Matrix{Sad.FloatM}, H_sparse),
                convert(Matrix{Sad.FloatM}, W_sparse),
                copy(S),
            )
            @test size(H_out, 1) < size(H_sparse, 1)
            @test size(H_out, 1) == size(W_out, 1)
            @test x_out[1] == 0.0
        end
    end
end

# B-spline and preprocess tests

@testset "B-spline smoothing" begin
    @testset "spline_design_matrix" begin
        x = collect(0.0:0.5:10.0)
        B = BSplineKit.BSplineBasis(4, collect(range(0, 10, length=8)))
        A = Sad.spline_design_matrix(x, B)
        @test size(A) == (length(x), length(B))
    end

    @testset "penalty_matrix" begin
        D = Sad.penalty_matrix(10, 2)
        @test size(D) == (8, 10)
        @test norm(D * ones(10)) < 1e-10
    end

    @testset "gcv_select" begin
        x = collect(0.0:0.5:10.0)
        y = 1.0 .+ 2.0 .* x .+ 0.5 .* sin.(x)
        B = BSplineKit.BSplineBasis(4, collect(range(0, 10, length=8)))
        lam = Sad.gcv_select(x, y, B)
        @test lam > 0
        @test isfinite(lam)
    end

    @testset "smooth_profile_bspline" begin
        x = collect(0.0:500.0:10000.0)
        y_true = 100.0 .+ 5e-5 .* x
        y_noisy = y_true .+ 0.1 .* randn(length(x))

        spl = Sad.smooth_profile_bspline(x, y_noisy; nknots=10, lambda=1.0)
        y_smooth = spl.(x)
        @test length(y_smooth) == length(x)
        @test maximum(abs.(y_smooth .- y_true)) < 2.0

        spl_gcv = Sad.smooth_profile_bspline(x, y_noisy; nknots=10)
        y_gcv = spl_gcv.(x)
        @test length(y_gcv) == length(x)
    end

    @testset "slope_from_spline" begin
        x = collect(0.0:500.0:10000.0)
        y_true = 100.0 .+ 5e-5 .* x
        spl = Sad.smooth_profile_bspline(x, y_true; nknots=10, lambda=0.01)
        slopes = Sad.slope_from_spline(spl, x)
        @test all(slopes .> 0)
        @test abs(mean(slopes) - 5e-5) < 2e-4
    end

    @testset "enforce_monotonicity" begin
        x = collect(0.0:500.0:10000.0)
        y = 100.0 .+ 5e-5 .* x .+ 0.5 .* sin.(2*pi .* x ./ 5000.0)
        spl = Sad.smooth_profile_bspline(x, y; nknots=10, lambda=0.01)
        spl_mono = Sad.enforce_monotonicity(spl, x, 1e-5)
        slopes = Sad.slope_from_spline(spl_mono, x)
        @test all(slopes .>= 1e-5)
    end
end

@testset "preprocess" begin
    reach, p = begin
        Dataset(TESTDATA) do f
            g = f.group["XS_Timeseries"]
            ri = f.group["River_Info"]
            x = (g["X"][:][end] .- g["X"][:])[end:-1:1, 1]
            H = convert(Matrix{Sad.FloatM}, g["H"][end:-1:1, 1:10:end])
            W = convert(Matrix{Sad.FloatM}, g["W"][end:-1:1, 1:10:end])
            S = (diff(H, dims=1) ./ diff(x))[:, 1:10:end]
            S_row = S[1:1, :]
            S = convert(Matrix{Sad.FloatM}, vcat(S_row, S))
            reach = Sad.preprocess(x, H, W, S)
            p = Sad.priors(Float64(ri["QWBM"][1]), reach.hmin, Sad.sinuous; reach)
            reach, p
        end
    end

    Smat = begin
        Dataset(TESTDATA) do f
            g = f.group["XS_Timeseries"]
            x = (g["X"][:][end] .- g["X"][:])[end:-1:1, 1]
            H = convert(Matrix{Sad.FloatM}, g["H"][end:-1:1, 1:10:end])
            S = (diff(H, dims=1) ./ diff(x))[:, 1:10:end]
            S_row = S[1:1, :]
            convert(Matrix{Sad.FloatM}, vcat(S_row, S))
        end
    end

    @testset "preprocess creates valid reach" begin
        reach_bsp = Sad.preprocess(reach.obs.x, reach.obs.H, reach.obs.W, Smat;
                                      dx=200.0, min_slope=1e-5)
        @test reach_bsp.nx == reach.nx
        @test reach_bsp.nt == reach.nt
        @test reach_bsp.nobs == reach.nobs
        @test length(reach_bsp.x) == reach_bsp.nx
        @test size(reach_bsp.H) == (reach_bsp.nx, reach_bsp.nt)
        @test size(reach_bsp.W) == (reach_bsp.nx, reach_bsp.nt)
        # S0 should match v1 (same method)
        @test abs(reach_bsp.S0(reach_bsp.x[1]) - reach.S0(reach.x[1])) < 1e-8
        @test abs(reach_bsp.S0(reach_bsp.x[end]) - reach.S0(reach.x[end])) < 1e-8
        # Valid timesteps should match
        @test length(findall(reach_bsp.valid)) == length(findall(reach.valid))
    end

    @testset "preprocess without Sobs" begin
        reach_bsp_noS = Sad.preprocess(reach.obs.x, reach.obs.H, reach.obs.W;
                                          dx=200.0, min_slope=1e-5)
        @test reach_bsp_noS.nx > 0
        @test reach_bsp_noS.nt == reach.nt
        @test minimum(reach_bsp_noS.S0.(reach_bsp_noS.x)) >= 1e-5
    end

    @testset "preprocess + infer" begin
        reach_bsp = Sad.preprocess(reach.obs.x, reach.obs.H, reach.obs.W, Smat;
                                      dx=200.0, min_slope=1e-5)
        res_bsp = Sad.infer(p, reach_bsp; σ_obs=1.0, ν=3.0, λ_smooth=1.0, iterations=100)
        @test length(res_bsp.Q_post) == reach_bsp.nt
        @test all(res_bsp.Q_post[findall(reach_bsp.valid)] .> 0)
        @test res_bsp.n_post > 0
        @test res_bsp.r_post > 0
    end
end