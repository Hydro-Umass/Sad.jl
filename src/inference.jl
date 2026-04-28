using LinearAlgebra
using Statistics
using StatsBase
using EmpiricalDistributions
using ProgressMeter

"""
    compute_A0(reach, reach_ensemble) -> Float64

Compute the reference cross-sectional area A0 at the minimum observed
downstream WSE (reach.hmin) using the posterior mean Dingman parameters.

In the Durand-Manning formulation:
    Q = (1/n) * (A0 + dA)^(5/3) * W^(-2/3) * S^(1/2)

A0 is the cross-sectional area at the minimum observed WSE — the
reference state from which SWOT dA observations are measured.

Uses the Dingman power-law cross section:
    A = Wb * (Ym / Yb)^(1/r) * y

where y is the mean flow depth at hmin and Ym = (r+1)/r * y.
"""
function compute_A0(reach::SWOTReach, reach_ensemble::Matrix{Float64})
    r_mean  = mean(reach_ensemble[2, :])
    z0_mean = mean(reach_ensemble[3, :])
    α_mean  = mean(reach_ensemble[4, :])

    # downstream node geometry
    x_ds = reach.x[1]
    Wb   = reach.wbf(x_ds)                      # bankfull width [m]
    hbf  = reach.hbf(x_ds)                      # bankfull WSE [m]
    zx   = z0_mean + reach.z(x_ds) * α_mean     # bed elevation at x_ds
                                                  # = z0_mean since z(x[1])=0
    Yb   = max(hbf - zx, 0.01)                  # bankfull mean depth [m]

    # mean flow depth at minimum observed WSE
    y0   = max((reach.hmin - zx) * r_mean / (r_mean + 1), 0.01)

    # Dingman area at y0
    Ym   = (r_mean + 1) / r_mean * y0
    A0   = Wb * (Ym / Yb)^(1 / r_mean) * y0

    return A0
end

"""
    letkf(A, d, HA, R; xs, ys, ρ) -> Matrix{Float64}

Local Ensemble Transform Kalman Filter.

State vector rows: [log(Q), n, r, z0, α]
Q is assimilated in log-space to enforce Q > 0 after the update and
improve Gaussianity of the Q marginal distribution.

# Arguments
- `A`:   (ndim × nens) state ensemble — row 1 is log(Q)
- `d`:   (nobs,) observation vector
- `HA`:  (nobs × nens) model-predicted observation ensemble
- `R`:   (nobs,) diagonal observation error variances
- `xs`:  state index patches for localisation (default: global)
- `ys`:  observation index patches for localisation (default: global)
- `ρ`:   covariance inflation factor (default 1.05)
"""
function letkf(A::Matrix{Float64}, d::Vector{Float64},
               HA::Matrix{Float64}, R::Vector{Float64};
               xs::Union{Nothing, Vector{Vector{Int}}} = nothing,
               ys::Union{Nothing, Vector{Vector{Int}}} = nothing,
               ρ::Float64 = 1.05)
    ndim, nens = size(A)
    nobs       = length(d)
    Aa         = zeros(ndim, nens)

    R_mat = Diagonal(R)

    Y = HA .- mean(HA, dims=2)
    X = A  .- mean(A,  dims=2)

    if isnothing(xs) || isnothing(ys)
        xs = [collect(1:ndim)]
        ys = [collect(1:nobs)]
    end

    for (lx, ly) in zip(xs, ys)
        Xl = X[lx, :]
        Yl = Y[ly, :]
        Rl = R_mat[ly, ly]

        C  = Yl' * inv(Matrix(Rl))
        P  = inv((nens - 1) / ρ * I + C * Yl)
        W  = real(sqrt(Symmetric((nens - 1) * P)))
        w  = P * C * (d[ly] .- mean(HA[ly, :], dims=2)[:])
        W  = W .+ reshape(w, :, 1)

        Aa[lx, :] = Xl * W .+ mean(A[lx, :], dims=2)
    end

    return Aa
end

"""
    infer(priors, reach; kwargs...) -> NamedTuple

Full SAD inference pipeline for a single SWOT reach.

Two-stage approach:
1. `infer_channel_params`: SIES with mean Q to estimate static hydraulic parameters [n, r, z0, α].
2. `infer_discharge`: LETKF per timestep to estimate dynamic Q conditioned on the parameter posterior.

# Arguments
- `p`:        `SWOTPriors` with prior distributions for all parameters.
- `reach`:    `SWOTReach` with observations and reach geometry.
- `time_str`: Optional timestamps used to select monthly discharge priors.
- `N`:        Ensemble size (used in both stages).
- `σₒ`:       Observation noise std [m].
- `γ`:        SIES damping factor (stage 1).
- `max_iter`: SIES maximum iterations (stage 1).
- `tol`:      SIES convergence tolerance (stage 1).
- `rho`:      LETKF covariance inflation factor (stage 2).

# Returns
A `NamedTuple` with:
- `Q_post`:         `Vector{Float64}` — posterior mean Q per timestep (NaN if invalid).
- `reach_ensemble`: `Matrix{Float64}` — 4×N posterior parameter ensemble [n, r, z0, α].
- `A_post`:         per-timestep state ensemble from LETKF (or `nothing` for invalid timesteps).
- `valid_ts`:       `Vector{Int}` — timestep indices with ≥ 2 valid observations.
- `completeness`:   `nothing` (retained for backwards compatibility).
"""
function infer(p::SWOTPriors, reach::SWOTReach;
               time_str::Vector{String} = String[],
               N::Int        = 1000,
               σₒ::Float64   = 0.5,
               rho::Float64  = 1.05)
    valid_ts = findall(reach.valid)
    if isempty(valid_ts)
        @warn "infer: no valid timesteps — returning empty posterior"
        return (reach_ensemble = Matrix{Float64}(undef, 4, 0),
                Q_post         = fill(NaN, reach.nt),
                A_post         = nothing,
                valid_ts       = valid_ts,
                completeness   = nothing)
    end

    # Stage 1: static channel parameters via SIES
    months = if !isempty(time_str) && p.Qp isa Vector
        map(time_str) do s
            try Month(DateTime(s)).value
            catch; 1
            end
        end
    else
        fill(1, reach.nt)
    end
    pa = infer_channel_params(p, reach, months,  N, σₒ)

    # Stage 2: dynamic discharge via LETKF
    result = infer_discharge(pa, reach;
                             sigma_obs=σₒ, N=N, time_str=time_str, rho=rho)
    reach_ens = hcat(rand(pa.np, N), rand(pa.rp, N), rand(pa.zp, N), rand(pa.ap, N))' |> Matrix
    return (reach_ensemble = reach_ens,
            Q_post         = result.Q_post,
            A_post         = result.A_post,
            valid_ts       = valid_ts,
            completeness   = nothing)
end

function infer_channel_params(p::SWOTPriors, reach::SWOTReach, months, N::Int = 1000, σₒ::Float64 = 0.5)
    logit(u) = log(u / (1 - u))
    logit1(u) = 1 / (1 + exp(-u))
    zbnds = [minimum(p.zp), maximum(p.zp)]
    nbnds = [minimum(p.np), maximum(p.np)]
    r_lb = minimum(p.rp)
    X = zeros(4, N)
    X[1, :] = logit.((rand(p.np, N) .- nbnds[1]) ./ (nbnds[2] - nbnds[1]))
    X[2, :] = log.(rand(p.rp, N) .- r_lb)
    X[3, :] = logit.((rand(p.zp, N) .- zbnds[1]) ./ (zbnds[2] - zbnds[1]))
    X[4, :] = rand(p.ap, N)

    rep_ts = select_representative_timesteps(reach; n_bins=5)
    obs_per_t = map(rep_ts) do t
        good = findall(.!ismissing.(reach.obs.H[:, t]))
        (x = collect(reach.obs.x[good]),
         H = Float64.(reach.obs.H[good, t]),
         H_bc = reach.H[1, t],
         month = months[t])
    end
    d = vcat([o.H[2:end] for o in obs_per_t]...)
    HX = zeros(length(d), N)
    valid_ens = fill(true, N)
    for e in 1:N
        row = 1
        for o in obs_per_t
            Qp_t = monthly_q_prior(p, o.month)
            Q_t = mean(Qp_t.untruncated)
            n_e = nbnds[1] + (nbnds[2] - nbnds[1]) * logit1(X[1, e])
            r_e = exp(X[2, e]) + r_lb
            z0_e = zbnds[1] + (zbnds[2] - zbnds[1]) * logit1(X[3, e])
            a_e = X[4, e]
            pred = gvf_solve(Q_t, n_e, r_e, a_e, z0_e, o.H_bc, reach;
                             saveat=o.x)
            if isnothing(pred) || !all(isfinite.(pred))
                HX[row:row+length(o.H)-2, e] .= NaN
                valid_ens[e] = false
            else
                HX[row:row+length(o.H)-2, e] = pred[2:end]
            end
            row += length(o.H) - 1
        end
    end

    # d = mean.(skipmissing.(eachrow(reach.obs.H)))
    valid_obs = any.(.!ismissing.(eachcol(reach.obs.H))) |> sum
    # nobs = length(d)
    # qm = p.Qp isa Vector ? mean([mean(qp.untruncated) for qp in p.Qp]) : mean(p.Qp.untruncated)
    # HX = zeros(nobs-1, N)
    # valid_ens = fill(true, N)
    # for e in 1:N
    #     n_e = X[1, e]
    #     r_e = X[2, e]
    #     z0_e = zbnds[1] + (zbnds[2] - zbnds[1]) * logit1(X[3, e])
    #     a_e = X[4, e]
    #     pred = gvf_solve(qm, n_e, r_e, a_e, z0_e, d[1], reach; saveat=reach.obs.x)
    #     if isnothing(pred) || !all(isfinite.(pred))
    #         valid_ens[e] = false
    #         continue
    #     end
    #     HX[:, e] = pred[2:end]
    # end
    good_idx = findall(valid_ens)
    HX = HX[:, good_idx]
    X = X[:, good_idx]

    R = fill(σₒ / sqrt(valid_obs), length(d[2:end])).^2
    Xa = letkf(X, d[2:end], HX, R)
    return SWOTPriors(p.Qp,
                      UvBinnedDist(fit(Histogram, nbnds[1] .+ (nbnds[2] - nbnds[1]) .* logit1.(Xa[1, :]))),
                      UvBinnedDist(fit(Histogram, exp.(Xa[2, :]) .+ r_lb)),
                      UvBinnedDist(fit(Histogram, zbnds[1] .+ (zbnds[2] - zbnds[1]) .* logit1.(Xa[3, :]))),
                      UvBinnedDist(fit(Histogram, Xa[4, :])))
end

"""
    infer_discharge(p, reach; kwargs...) -> NamedTuple

Stage-2 LETKF discharge inference. For each valid timestep, samples hydraulic
parameters from `p`'s parameter distributions (typically posterior empirical
distributions from `infer_channel_params`) and draws Q from the (monthly) prior,
then runs the LETKF update against observed WSE.

# Arguments
- `p`:     `SWOTPriors` — parameter distributions (may be posterior from stage 1).
- `reach`: `SWOTReach`.

# Keyword arguments
- `sigma_obs`: Observation noise std [m] (default 0.1).
- `N`:         Ensemble size per timestep (default 500).
- `time_str`:  Timestamps for monthly prior selection.
- `rho`:       LETKF covariance inflation factor (default 1.05).

# Returns
`NamedTuple` with `Q_post` (length nt) and `A_post` (per-timestep state ensemble).
"""
function infer_discharge(p::SWOTPriors, reach::SWOTReach;
                         sigma_obs::Float64       = 0.1,
                         N::Int                   = 500,
                         time_str::Vector{String} = String[],
                         rho::Float64             = 1.05)
    months = if !isempty(time_str) && p.Qp isa Vector
        map(time_str) do s
            try Month(DateTime(s)).value
            catch; 1
            end
        end
    else
        fill(1, reach.nt)
    end
    if sum(reach.valid) == 0
        @warn "No valid timesteps: returning empty posterior"
        return (Q_post = fill(NaN, reach.nt),
                A_post = Vector{Union{Nothing, Matrix{Float64}}}(nothing, reach.nt))
    end
    nt     = reach.nt
    Q_post = fill(NaN, nt)
    A_post = Vector{Union{Nothing, Matrix{Float64}}}(nothing, nt)
    valid_ts = findall(reach.valid)
    n_est = mean(p.np)
    r_est = mean(p.rp)
    z0_est = mean(p.zp)
    a_est = mean(p.ap)
    n_ens = rand(p.np, N)
    r_ens = rand(p.rp, N)
    z0_ens = rand(p.zp, N)
    a_ens = rand(p.ap, N)
    @showprogress "Timesteps" for t in valid_ts
        Q_ens  = rand(monthly_q_prior(p, months[t]), N)
        # A    = reshape(Q_ens, 1, N)
        A    = reshape(log.(Q_ens), 1, N)
        H_bc = reach.H[1, t]
        good = findall(.!ismissing.(reach.obs.H[:, t]))
        if length(good) < 2
            @warn "Skipping timestep $t: fewer than 2 observations"
            continue
        end
        d = Float64.(reach.obs.H[good, t])
        x_obs = collect(reach.obs.x[good])
        n_q   = length(Q_ens)
        HA = zeros(length(d), n_q)
        good_ens = Int[]
        for k in 1:n_q
            # pred = gvf_solve((A[1, k]), n_est, r_est, a_est, z0_est, H_bc, reach; saveat=x_obs)
            pred = gvf_solve(exp(A[1, k]), n_ens[k], r_ens[k], a_ens[k], z0_ens[k], H_bc, reach; saveat=x_obs)
            # pred = gvf_solve((A[1, k]), n_ens[k], r_ens[k], a_ens[k], z0_ens[k], H_bc, reach; saveat=x_obs)
            if !isnothing(pred) && all(isfinite.(pred))
                HA[:, k] = pred
                push!(good_ens, k)
            end
        end
        if length(good_ens) < 2
            @warn "Skipping timestep $t: fewer than 2 valid ensemble members"
            continue
        end
        A_filt  = A[:, good_ens]
        HA_filt = HA[:, good_ens]
        sigma_node = fill(sigma_obs, length(d)-1)
        R_diag = sigma_node .^ 2
        A_analysis = letkf(A_filt, d[2:end], HA_filt[2:end, :], R_diag; ρ=rho)
        Qp_t = monthly_q_prior(p, months[t])
        Q_lo = quantile(Qp_t, 0.001)
        Q_hi = quantile(Qp_t, 0.999)
        Q_post_ens = clamp.(exp.(A_analysis[1, :]), Q_lo, Q_hi)
        # Q_post_ens = clamp.(A_analysis[1, :], Q_lo, Q_hi)
        if !all(isfinite.(Q_post_ens))
            @warn "Timestep $t: non-finite Q after exponentiation — skipping"
            continue
        end
        Q_post[t] = median(Q_post_ens)
        A_natural      = copy(A_analysis)
        A_natural[1,:] = Q_post_ens
        A_post[t]      = A_natural
    end
    return (
        Q_post         = Q_post,
        A_post = A_post,
    )
end
