using LinearAlgebra
using Statistics
using ProgressMeter

"""
    obs_completeness(reach) -> Vector{Float64}

Compute temporal data completeness for each SWOT node — the fraction of
timesteps with a valid (non-missing) H observation. Used to inflate
observation error for nodes that are frequently missing.
"""
function obs_completeness(reach::SWOTReach)
    nobs, nt = size(reach.obs.H)
    valid_frac = [count(.!ismissing.(reach.obs.H[j, :])) / nt for j in 1:nobs]
    return clamp.(valid_frac, 1.0 / nt, 1.0)
end

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

    # Downstream node geometry
    x_ds = reach.x[1]
    Wb   = reach.wbf(x_ds)                      # bankfull width [m]
    hbf  = reach.hbf(x_ds)                      # bankfull WSE [m]
    zx   = z0_mean + reach.z(x_ds) * α_mean     # bed elevation at x_ds
                                                  # = z0_mean since z(x[1])=0
    Yb   = max(hbf - zx, 0.01)                  # bankfull mean depth [m]

    # Mean flow depth at minimum observed WSE
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

# Arguments
- `priors`:               SWOTPriors
- `reach`:                SWOTReach
- `sigma_obs`:            observation error std dev [m] (default 0.1)
- `N`:                    rejection sample target size (default 500)
- `eps_rel`:              rejection RMSE threshold (default 0.5)
- `eps_abs_q`:            upstream WSE tolerance for q_ensemble [m] (default 1.0)
- `n_bins`:               quantile bins for representative timesteps (default 5)
- `rho`:                  LETKF covariance inflation (default 1.05)
- `completeness_weight`:  inflate R by inverse completeness (default true)

# Returns
NamedTuple with fields:
- `reach_ensemble`: (4 × N_acc) matrix [n, r, z0, α]
- `Q_post`:         (nt,) vector of posterior mean Q (NaN where skipped)
- `A_post`:         vector of posterior ensemble matrices (nothing where skipped)
- `rep_ts`:         representative timestep indices
- `completeness`:   (nobs,) per-node data completeness
"""
function infer(priors::SWOTPriors, reach::SWOTReach;
               sigma_obs::Float64           = 0.1,
               N::Int                       = 500,
               eps_rel::Float64             = 0.5,
               eps_abs_q::Float64           = 1.0,
               n_bins::Int                  = 5,
               rho::Float64                 = 1.05,
               completeness_weight::Bool    = true,
               time_str::Vector{String}     = String[])

    # parse months from time strings if provided and Qp is monthly
    months = if !isempty(time_str) && priors.Qp isa Vector
        map(time_str) do s
            try Month(DateTime(s)).value
            catch; 1
            end
        end
    else
        fill(1, reach.nt)
    end

    # refine reach parameter prior
    # return early if no valid timesteps, e.g. all observations missing
    if sum(reach.valid) == 0
        @warn "No valid timesteps: returning empty posterior"
        return (
            reach_ensemble = Matrix{Float64}(undef, 4, 0),
            Q_post         = fill(NaN, reach.nt),
            A_post         = Vector{Union{Nothing, Matrix{Float64}}}(nothing, reach.nt),
            rep_ts         = Int[],
            completeness   = Float64[],
        )
    end
    rep_ts         = select_representative_timesteps(reach; n_bins=n_bins)
    reach_ensemble = rejection_sample(priors, reach, months;
                                      N=N, ε_rel=eps_rel, n_bins=n_bins)
    if size(reach_ensemble, 2) == 0
        error("Rejection sampling produced no accepted samples.")
    end

    N_acc = size(reach_ensemble, 2)
    @info "Reach ensemble size: $N_acc"

    completeness = obs_completeness(reach)
    @info "Node completeness: min=$(round(minimum(completeness), digits=2))  " *
          "median=$(round(median(completeness), digits=2))  " *
          "max=$(round(maximum(completeness), digits=2))"

    nt     = reach.nt
    Q_post = fill(NaN, nt)
    A_post = Vector{Union{Nothing, Matrix{Float64}}}(nothing, nt)

    valid_ts = findall(reach.valid)

    @showprogress "Timesteps" for t in valid_ts

        # --- Stage 2: Q ensemble conditioned on upstream WSE ---
        Q_ens, member_idx = q_ensemble(priors, reach_ensemble, reach, t,
                                        months[t];
                                        eps_abs=eps_abs_q)
        if isempty(Q_ens)
            @warn "Skipping timestep $t: empty Q ensemble"
            continue
        end

        n_ens  = reach_ensemble[1, member_idx]
        r_ens  = reach_ensemble[2, member_idx]
        z0_ens = reach_ensemble[3, member_idx]
        a_ens  = reach_ensemble[4, member_idx]

        # State ensemble in log(Q) space — rows: [log(Q), n, r, z0, α]
        A    = vcat(log.(Q_ens)', n_ens', r_ens', z0_ens', a_ens')
        H_bc = reach.H[1, t]

        # --- Stage 3: forward model ---
        good = findall(.!ismissing.(reach.obs.H[:, t]))
        if length(good) < 2
            @warn "Skipping timestep $t: fewer than 2 observations"
            continue
        end

        d     = Float64.(reach.obs.H[good, t])
        x_obs = collect(reach.obs.x[good])
        n_q   = length(Q_ens)

        HA = zeros(length(d), n_q)
        good_ens = Int[]
        for k in 1:n_q
            pred = gvf_solve(exp(A[1, k]), A[2, k], A[3, k], A[5, k], A[4, k],
                             H_bc, reach; saveat=x_obs)
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

        # --- Stage 4: build R ---
        if completeness_weight
            sigma_node = sigma_obs ./ sqrt.(completeness[good])
        else
            sigma_node = fill(sigma_obs, length(d))
        end
        R_diag = sigma_node .^ 2

        # --- Stage 5: LETKF update in log(Q) space ---
        A_analysis = letkf(A_filt, d, HA_filt, R_diag; ρ=rho)

        # Convert log(Q) back to natural space
        Q_post_ens = exp.(A_analysis[1, :])

        if !all(isfinite.(Q_post_ens))
            @warn "Timestep $t: non-finite Q after exponentiation — skipping"
            continue
        end

        Q_post[t] = mean(Q_post_ens)

        # Store posterior in natural space [Q, n, r, z0, α]
        A_natural      = copy(A_analysis)
        A_natural[1,:] = Q_post_ens
        A_post[t]      = A_natural

        @info "Timestep $t: Q = $(round(Q_post[t], digits=1)) m³/s  " *
              "spread = $(round(std(Q_post_ens), digits=1))"
    end

    return (
        reach_ensemble = reach_ensemble,
        Q_post         = Q_post,
        A_post         = A_post,
        rep_ts         = rep_ts,
        completeness   = completeness,
    )
end
