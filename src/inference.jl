using LinearAlgebra
using Statistics
using Optim
using ForwardDiff

# ============================================================
# Joint MAP estimation with Manning + GVF perturbation
# ============================================================

"""
    safe_logpdf(d, x)

Compute log-pdf with smooth barrier for bounded distributions.

For distributions with bounded support (e.g., Uniform), `-Inf` values outside
the support crash ForwardDiff-based optimizers. This function replaces bounded
log-pdf with a smooth approximation: a Normal matching the distribution's mean
and variance, which has infinite support and is always ForwardDiff-differentiable.

For distributions with unbounded support (Normal, LogNormal, etc.), `logpdf`
is returned unchanged.
"""
function safe_logpdf(d, x)
    if d isa Uniform
        # Replace Uniform(μ, σ) → Normal(mean, std) matching first two moments
        μ_u = (minimum(d) + maximum(d)) / 2
        σ_u = (maximum(d) - minimum(d)) / (2 * sqrt(3))
        return logpdf(Normal(μ_u, σ_u), x)
    elseif d isa Truncated
        # Use the untruncated distribution logpdf but add a smooth quadratic
        # barrier outside the truncation bounds. This preserves the
        # distribution shape within bounds while providing a strong smooth
        # gradient back toward the support when the optimizer strays outside.
        lp = safe_logpdf(d.untruncated, x)
        lo, hi = minimum(d), maximum(d)
        # Barrier width: fraction of the support range
        barrier_scale = (hi - lo) * 0.2
        if x < lo
            lp -= 0.5 * ((lo - x) / barrier_scale)^2
        elseif x > hi
            lp -= 0.5 * ((x - hi) / barrier_scale)^2
        end
        return lp
    else
        return logpdf(d, x)
    end
end

"""
    ManningPrecomp

Pre-computed reach data for the Manning forward model. Stores
geometry and observation data that are constant across optimization
iterations, avoiding redundant recomputation.

# Fields
- `x_nodes`:      SWOT observation node chainage [m]
- `Wb_nodes`:     bankfull width at each observation node [m]
- `hbf_z_nodes`:  hbf(x) − z_ref(x) at each node [m] — the z₀-independent
                   part of bankfull depth; Yb(z₀) = hbf_z_nodes − z₀
- `z_ref_nodes`:  cumulative bed elevation reference z_ref at each node [m]
              (z_ref[1] = 0 by convention; z₀ + z_ref(x) = total bed elevation)
- `S0`:           reach-averaged bed slope [-]
- `nt`:           total number of timesteps
- `valid_ts`:     indices of valid timesteps
- `H_bc`:         downstream boundary WSE, one per valid timestep [m]
- `upstream_j`:   for each valid timestep, 1-based indices (into `x_nodes`)
                   of non-missing upstream observations (j ≥ 2)
- `upstream_H`:   for each valid timestep, observed WSE at those nodes [m]
- `months_v`:     month (1–12) for each valid timestep
"""
struct ManningPrecomp
    x_nodes::Vector{Float64}
    Wb_nodes::Vector{Float64}
    hbf_z_nodes::Vector{Float64}
    z_ref_nodes::Vector{Float64}
    S0::Float64
    nt::Int
    valid_ts::Vector{Int}
    H_bc::Vector{Float64}
    upstream_j::Vector{Vector{Int}}
    upstream_H::Vector{Vector{Float64}}
    months_v::Vector{Int}
    σ_obs_est::Float64  # automatically estimated WSE observation noise [m]
    upstream_W::Vector{Vector{Float64}}  # observed width at upstream nodes
    σ_W_est::Float64    # automatically estimated width observation noise [m]
end

"""
    precompute_manning(reach, months; S0) -> ManningPrecomp

Pre-compute reach geometry and observation data for the Manning
forward model. Call once before optimization.

# Arguments
- `reach`:  `SWOTReach` from preprocessing
- `months`: vector of month indices (1–12), one per timestep
- `S0`:     optional reach-averaged bed slope; if `nothing`, computed as
            `mean(reach.S0.(reach.x))`
"""
function precompute_manning(reach::SWOTReach, months::Vector{Int};
                            S0::Union{Nothing, Float64} = nothing)
    x_nodes    = collect(reach.obs.x)
    Wb_nodes   = reach.wbf.(x_nodes)
    hbf_nodes  = reach.hbf.(x_nodes)
    z_ref_nodes = reach.z.(x_nodes)
    hbf_z_nodes = hbf_nodes .- z_ref_nodes   # constant part of Yb(z₀)
    S0_val = isnothing(S0) ? mean(reach.S0.(reach.x)) : S0

    valid_ts = findall(reach.valid)
    nt = reach.nt

    H_bc_list    = Float64[]
    upstream_j   = Vector{Int}[]
    upstream_H   = Vector{Float64}[]
    upstream_W   = Vector{Float64}[]
    months_v     = Int[]

    for t in valid_ts
        # Downstream boundary condition: use the most downstream RAW observation
        # rather than the interpolated grid, which can be unreliable when the
        # downstream node has missing data and PCHIP extrapolates poorly.
        H_bc_val = NaN
        for j in 1:length(x_nodes)
            if !ismissing(reach.obs.H[j, t])
                H_bc_val = Float64(reach.obs.H[j, t])
                break
            end
        end
        # Fallback to interpolated grid if no raw observation exists
        if isnan(H_bc_val)
            H_bc_val = reach.H[1, t]
        end
        push!(H_bc_list, H_bc_val)
        # Collect non-missing upstream observations (j >= 2)
        js   = Int[]
    Hobs = Float64[]
    Wobs = Float64[]
        for j in 2:length(x_nodes)
            h_valid = !ismissing(reach.obs.H[j, t])
            w_valid = !ismissing(reach.obs.W[j, t])
            if h_valid || w_valid
                push!(js, j)
                push!(Hobs, h_valid ? Float64(reach.obs.H[j, t]) : NaN)
                push!(Wobs, w_valid ? Float64(reach.obs.W[j, t]) : NaN)
            end
        end
        push!(upstream_j, js)
        push!(upstream_H, Hobs)
        push!(upstream_W, Wobs)
        push!(months_v, months[t])
    end

    # Estimate WSE observation noise from the data
    σ_obs_est = estimate_σ_obs(x_nodes, upstream_j, upstream_H, S0_val)

    # Estimate width observation noise from the data
    σ_W_est = estimate_σ_W(upstream_W)

    ManningPrecomp(x_nodes, Wb_nodes, hbf_z_nodes, z_ref_nodes, S0_val, nt,
                   valid_ts, H_bc_list, upstream_j, upstream_H, months_v, σ_obs_est,
                   upstream_W, σ_W_est)
end

"""
    estimate_σ_obs(x_nodes, upstream_j, upstream_H, S0) -> Float64

Estimate WSE observation noise from the raw SWOT observations.

For each valid timestep, compute the standard deviation of residuals
between observed WSE and a smooth monotonic profile (approximated by the
expected WSE from bed slope). The median across timesteps gives a robust
noise estimate. A minimum floor of 0.3 m is applied to account for model
structural error beyond pure observation noise.

# Arguments
- `x_nodes`:   observation node chainage [m]
- `upstream_j`: per-timestep indices of non-missing upstream nodes
- `upstream_H`: per-timestep observed WSE at those nodes [m]
- `S0`:         reach-averaged bed slope [-]

# Returns
Estimated WSE observation error σ_obs [m]
"""
function estimate_σ_obs(x_nodes::Vector{Float64},
                        upstream_j::Vector{Vector{Int}},
                        upstream_H::Vector{Vector{Float64}},
                        S0::Float64)
    noise_ests = Float64[]
    for k in 1:length(upstream_j)
        js = upstream_j[k]
        Hobs = upstream_H[k]
        x_vals = x_nodes[js]
        # Filter to only non-NaN H obs (node may have width but not WSE)
        valid_idx = findall(idx -> !isnan(Hobs[idx]), 1:length(js))
        if length(valid_idx) < 3
            continue
        end
        H_valid = Hobs[valid_idx]
        resid_sq = Float64[]
        for ii in 2:length(H_valid)-1
            r = H_valid[ii] - 0.5 * (H_valid[ii-1] + H_valid[ii+1])
            push!(resid_sq, r^2)
        end
        if !isempty(resid_sq)
            # σ from local residuals (factor 2/3 corrects for 3-point averaging)
            σ_t = sqrt(mean(resid_sq) * 2.0 / 3.0)
            push!(noise_ests, σ_t)
        end
    end
    if isempty(noise_ests)
        return 1.0  # conservative default
    end
    σ_est = median(noise_ests)
    # Apply minimum floor to account for model structural error
    return max(σ_est, 0.3)
end

"""
    estimate_σ_W(upstream_W) -> Float64

Estimate width observation noise from the raw SWOT observations.
Computes the median coefficient of variation (CV) of width across valid
upstream nodes for each timestep, then scales by median width.
A minimum floor ensures model structural error is accounted for.
"""
function estimate_σ_W(upstream_W::Vector{Vector{Float64}})
    cv_ests = Float64[]
    for k in 1:length(upstream_W)
        Wobs = filter(w -> !isnan(w) && w > 0, upstream_W[k])
        if length(Wobs) < 3
            continue
        end
        # CV of width across nodes estimates spatial variability,
        # which is a lower bound on observational noise
        cv = std(Wobs) / mean(Wobs)
        push!(cv_ests, cv)
    end
    if isempty(cv_ests)
        return 30.0  # conservative default [m]
    end
    median_cv = median(cv_ests)
    # Typical SWOT width uncertainty is 10-20% of width
    σ_estimate = max(median_cv * 0.5, 0.10)  # at least 10% CV
    # Convert CV to absolute noise: use mean width as reference
    # (will be scaled per-node in the likelihood using Wb)
    return σ_estimate  # return as fractional uncertainty (unitless)
end

"""
    neg_log_joint(θ, precomp, priors; σ_obs, ν, λ_smooth) -> Real

Negative log-joint for MAP estimation over all parameters simultaneously.

The parameter vector θ encodes:
- `θ[1:nt]` = log(Q_t)  — discharge at each timestep (log-space)
- `θ[nt+1]` = log(n)    — Manning roughness (log-space)
- `θ[nt+2]` = log(r)    — Dingman shape exponent (log-space)
- `θ[nt+3]` = z₀        — downstream bed elevation (natural space)

The objective comprises four terms:
1. **WSE likelihood**: Student-t (ν < ∞) or Gaussian (ν = Inf) per-node
   residuals between predicted and observed WSE at upstream nodes (j ≥ 2).
2. **Q priors**: monthly climatological log-normal priors in log-Q space.
3. **Static priors**: on n, r, z₀ from the SoS database.
4. **Temporal smoothness**: penalty on successive log-Q differences.

# Arguments
- `θ`:        parameter vector of length `nt + 3`
- `precomp`:  `ManningPrecomp` with pre-computed reach data
- `priors`:   `SWOTPriors` with prior distributions

# Keyword arguments
- `σ_obs`:     observation error scale [m] (default `NaN`: auto-estimated from data)
- `ν`:         Student-t degrees of freedom; `Inf` = Gaussian (default `5.0`)
- `λ_smooth`:  temporal smoothness weight on log-Q (default `0.0`, disabled)

# Returns
Scalar negative log-joint value. Minimized by `infer`.

# Conventions (see AGENTS.md)
- Q and n are estimated in log-space to enforce positivity and improve conditioning.
- z₀ is in natural space (can be negative).
- The downstream boundary WSE (H_bc, node j=1) is an INPUT to the forward model,
  not an observation; the likelihood applies at nodes j ≥ 2 only.
"""
function neg_log_joint(θ::AbstractVector, precomp::ManningPrecomp,
                       priors::SWOTPriors;
                       σ_obs::Float64 = NaN,
                       ν::Float64    = 5.0,
                       λ_smooth::Float64 = 0.0,
                       use_width::Bool = false)
    # Auto-select σ_obs from estimated noise if not provided
    σ = isnan(σ_obs) ? precomp.σ_obs_est : σ_obs
    nt = precomp.nt
    logQ = view(θ, 1:nt)
    logn = θ[nt + 1]
    logr = θ[nt + 2]
    z0   = θ[nt + 3]

    n = exp(logn)
    r = exp(logr)

    # Bankfull depth at each node (depends on z₀)
    Yb_nodes = max.(precomp.hbf_z_nodes .- z0, 0.01)

    nll = zero(z0)  # accumulator with correct Dual type

    # Pre-compute month index for each timestep (avoids repeated linear search)
    month_map = zeros(Int, nt)
    for (ki, t) in enumerate(precomp.valid_ts)
        month_map[t] = precomp.months_v[ki]
    end
    # Invalid timesteps default to month 1
    for t in 1:nt
        if month_map[t] == 0
            month_map[t] = 1
        end
    end

    # --- 1. WSE log-likelihood at upstream nodes ---
    for k in 1:length(precomp.valid_ts)
        Q_t = exp(logQ[precomp.valid_ts[k]])
        WSE_pred = manning_wse_backwater(
            Q_t, n, r, z0, precomp.S0, precomp.H_bc[k],
            precomp.x_nodes, precomp.Wb_nodes, Yb_nodes, precomp.z_ref_nodes)

        js   = precomp.upstream_j[k]
        Hobs = precomp.upstream_H[k]
        Wobs = precomp.upstream_W[k]
        for idx in 1:length(js)
            # WSE likelihood (skip if H is NaN — missing observation)
            if !isnan(Hobs[idx])
                residual_H = Hobs[idx] - WSE_pred[js[idx]]
                if isinf(ν)
                    nll += 0.5 * (residual_H / σ)^2
                else
                    nll += ((ν + 1) / 2) * log(1 + (residual_H / σ)^2 / ν)
                end
            end
            # Width likelihood (skip if W is NaN — missing observation)
            if use_width && !isnan(Wobs[idx])
                # Predicted width: W_pred = Wb * (y/Yb)^(1/r)
                y_at_node = max((WSE_pred[js[idx]] - z0 - precomp.z_ref_nodes[js[idx]]) * r / (r + 1), 0.01)
                Yb_at_node = max(precomp.hbf_z_nodes[js[idx]] - z0, 0.01)
                W_pred = precomp.Wb_nodes[js[idx]] * (y_at_node / Yb_at_node)^(1/r)
                # Use large width uncertainty to reflect model structural error
                σ_W = max(precomp.σ_W_est * precomp.Wb_nodes[js[idx]], 0.3 * precomp.Wb_nodes[js[idx]])
                residual_W = Wobs[idx] - W_pred
                nll += 0.5 * (residual_W / σ_W)^2
            end
        end
    end

    # --- 2. Q priors (in log-Q space) ---
    for t in 1:nt
        q_prior = monthly_q_prior(priors, month_map[t])
        d = q_prior isa Truncated ? q_prior.untruncated : q_prior
        if d isa LogNormal
            # logQ ~ Normal(d.μ, d.σ) in log-space — exact, always finite
            nll -= logpdf(Normal(d.μ, d.σ), logQ[t])
        else
            # Fallback: use safe_logpdf with Jacobian correction
            nll -= safe_logpdf(q_prior, exp(logQ[t])) + logQ[t]
        end
    end

    # --- 3. Static parameter priors ---
    # n prior (in natural n-space, with Jacobian for logn parameterization)
    nll -= safe_logpdf(priors.np, n) + logn

    # r prior (in natural r-space, with Jacobian for logr parameterization)
    nll -= safe_logpdf(priors.rp, r) + logr

    # z₀ prior (natural space, no Jacobian)
    nll -= safe_logpdf(priors.zp, z0)

    # --- 4. Temporal smoothness on log-Q ---
    for t in 2:nt
        nll += λ_smooth * (logQ[t] - logQ[t - 1])^2
    end

    return nll
end

"""
    infer(priors, reach; kwargs...) -> NamedTuple

SAD inference: joint MAP estimation over discharge and channel parameters
using the analytical Manning + first-order GVF perturbation forward model.

All parameters are estimated simultaneously by minimizing a single objective
(the negative log-joint) with L-BFGS, then Laplace approximation provides
posterior uncertainty.

# Parameter vector θ (length nt + 3)
- `θ[1:nt]`   = log(Q_t) — discharge at each timestep
- `θ[nt+1]`  = log(n)    — Manning roughness
- `θ[nt+2]`  = log(r)    — Dingman shape exponent
- `θ[nt+3]`  = z₀        — downstream bed elevation

# Arguments
- `priors`: `SWOTPriors` with prior distributions
- `reach`:  `SWOTReach` with observations and geometry

# Keyword arguments
- `time_str`:    timestamps for monthly Q-prior selection
- `σ_obs`:       observation error scale [m] (default `NaN`: auto-estimated from data)
- `ν`:           Student-t d.f.; `Inf` = Gaussian (default `5.0`)
- `λ_smooth`:    temporal smoothness weight (default `0.1`)
- `iterations`:  max L-BFGS iterations (default 500)
- `g_tol`:       gradient convergence tolerance (default 1e-6)
- `S0`:          optional reach-averaged slope override

# Returns
`NamedTuple` with:
- `Q_post`:  posterior mean Q per timestep (MAP estimate)
- `Q_std`:   posterior std of Q (Laplace delta method)
- `n_post`:  MAP Manning n
- `r_post`:  MAP Dingman r
- `z0_post`: MAP downstream bed elevation
- `θ_map`:   full MAP parameter vector
- `Σ`:       posterior covariance matrix (Laplace)
- `result`:  Optim.jl optimization result
- `valid_ts`: indices of valid timesteps
- `precomp`:  `ManningPrecomp` used in inference
"""
function infer(p::SWOTPriors, reach::SWOTReach;
               time_str::Vector{String} = String[],
               σ_obs::Float64    = NaN,
               ν::Float64        = 5.0,
               λ_smooth::Float64 = 0.1,
               iterations::Int   = 500,
               g_tol::Float64   = 1e-6,
               S0::Union{Nothing, Float64} = nothing,
               use_width::Bool  = false)
    # Extract months
    months = if !isempty(time_str) && p.Qp isa Vector
        map(time_str) do s
            try Month(DateTime(s)).value
            catch; 1
            end
        end
    else
        fill(1, reach.nt)
    end

    # Pre-compute reach geometry and observations
    precomp = precompute_manning(reach, months; S0=S0)

    if isempty(precomp.valid_ts)
        @warn "infer: no valid timesteps"
        nt = reach.nt
        return (Q_post  = fill(NaN, nt),
                Q_std    = fill(NaN, nt),
                n_post   = NaN, r_post = NaN, z0_post = NaN,
                θ_map    = fill(NaN, nt + 3),
                Σ        = zeros(nt + 3, nt + 3),
                result   = nothing,
                valid_ts = Int[], precomp = precomp)
    end

    # Use data-estimated σ_obs if not explicitly provided
    σ_use = isnan(σ_obs) ? precomp.σ_obs_est : σ_obs
    println("  σ_obs = $(round(σ_use, digits=3)) m")

    # Initialize from prior medians
    θ0 = initialize_theta(p, precomp)

    # Objective closure
    obj(θ) = neg_log_joint(θ, precomp, p; σ_obs=σ_use, ν=ν, λ_smooth=λ_smooth, use_width=use_width)

    # In-place gradient via ForwardDiff
    function g!(G, θ)
        G .= ForwardDiff.gradient(obj, θ)
    end

    # L-BFGS optimization
    result = optimize(obj, g!, θ0, LBFGS(),
                      Optim.Options(iterations=iterations,
                                     g_tol=g_tol,
                                     allow_f_increases=true,
                                     show_trace=false))

    θ_map = Optim.minimizer(result)

    # Extract physical parameters
    nt = reach.nt
    Q_post  = exp.(θ_map[1:nt])
    n_post  = exp(θ_map[nt + 1])
    r_post  = exp(θ_map[nt + 2])
    z0_post = θ_map[nt + 3]

    # Laplace uncertainty
    Σ = laplace_uncertainty(obj, θ_map)
    # Q_std via delta method: σ_Q ≈ Q · σ_{logQ}
    Q_std = Q_post .* sqrt.(diag(Σ)[1:nt])

    return (Q_post  = Q_post,
            Q_std   = Q_std,
            n_post  = n_post,
            r_post  = r_post,
            z0_post = z0_post,
            θ_map   = θ_map,
            Σ       = Σ,
            result  = result,
            valid_ts = precomp.valid_ts,
            precomp  = precomp)
end

"""
    initialize_theta(priors, precomp) -> Vector{Float64}

Construct initial parameter vector θ₀ from prior medians.

For Q, uses a physics-informed estimate: given prior-medians for n, r, z₀,
estimates mean depth from the first valid timestep's observations, then
computes Q via Manning's equation. Falls back to the prior median if the
estimate is infeasible.
"""
function initialize_theta(p::SWOTPriors, precomp::ManningPrecomp)
    nt = precomp.nt

    # Static parameters from prior medians
    n_init = quantile(p.np, 0.5)
    r_init = quantile(p.rp, 0.5)
    z0_init = quantile(p.zp, 0.5)

    # Initialize Q from monthly prior medians (more robust than depth-based
    # estimation, especially for noisy reaches where observed WSE gives
    # unreliable depth estimates)
    logQ_init = fill(0.0, nt)
    for t in 1:nt
        ki = findfirst(==(t), precomp.valid_ts)
        mo = isnothing(ki) ? 1 : precomp.months_v[ki]
        q_prior = monthly_q_prior(p, mo)
        q_median = quantile(q_prior, 0.5)
        logQ_init[t] = log(max(q_median, 0.01))
    end

    return vcat(logQ_init, log(n_init), log(r_init), z0_init)
end

"""
    laplace_uncertainty(obj, θ_map) -> Matrix

Compute the Laplace approximation to the posterior covariance by inverting
the Hessian of the negative log-joint at the MAP estimate.

Σ = H⁻¹ where H = ∇² neg_log_joint(θ_map)

The Hessian is computed via ForwardDiff. If the inverse is not positive
definite (due to numerical issues), a small diagonal regularization is added.

# Arguments
- `obj`:   scalar objective function f(θ) → neg_log_joint
- `θ_map`: MAP parameter vector

# Returns
Posterior covariance matrix Σ of size (length(θ_map) × length(θ_map)).
"""
function laplace_uncertainty(obj, θ_map::AbstractVector)
    H = ForwardDiff.hessian(obj, θ_map)
    # Symmetrize (ForwardDiff hessian should already be symmetric,
    # but machine precision may cause slight asymmetry)
    H = Symmetric(H)
    try
        Σ = inv(H)
        # Check for negative variances (indicates non-convexity at MAP)
        if any(diag(Σ) .< 0)
            @warn "laplace_uncertainty: negative diagonal detected, adding regularization"
            Σ = inv(H + 1e-4 * I)
        end
        return Matrix(Σ)
    catch e
        @warn "laplace_uncertainty: Hessian inversion failed" exception=e
        return zeros(length(θ_map), length(θ_map))
    end
end