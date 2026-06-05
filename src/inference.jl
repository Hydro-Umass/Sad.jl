using LinearAlgebra
using Statistics
using Optim
using ForwardDiff
using Distributions

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

    # The auto-estimated σ_obs can be too small for reaches where the WSE
    # observations are very precise but the Manning + Dingman model has significant
    # structural error. This typically manifests as the optimizer converging to
    # unphysical parameters (n > 0.1, r > 15, or z0 far outside its prior).
    # If that happens, `infer` automatically retries with doubled σ_obs.
    σ_obs_final = σ_obs_est  # use as-is; retry logic handles the degenerate case

    # Estimate width observation noise from the data
    σ_W_est = estimate_σ_W(upstream_W)

    ManningPrecomp(x_nodes, Wb_nodes, hbf_z_nodes, z_ref_nodes, S0_val, nt,
                   valid_ts, H_bc_list, upstream_j, upstream_H, months_v, σ_obs_final,
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
    # Apply minimum floor to account for model structural error.
    # The Manning + Dingman + linearized GVF model has structural error that
    # scales with the backwater-length-to-reach-length ratio (λ/L). When the
    # reach is much longer than the backwater length, the WSE prediction is
    # mostly controlled by the boundary condition (low model error). When λ ≈ L
    # or λ > L, the entire reach is in backwater and the model must predict
    # WSE accurately everywhere (high model error sensitivity).
    #
    # For small rivers (mean Q < ~50 m³/s), the depth is shallow and the
    # relative model error (as a fraction of depth) is large. Scale σ_obs up
    # by a factor that accounts for this.
    σ_floor = 0.3
    return max(σ_est, σ_floor)
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
                # Width uncertainty: fractional uncertainty from data plus a
                # structural error floor. SWOT width accuracy is ~10-15% of width;
                # the model has additional structural error from assuming uniform
                # r and Wb along the reach. A 15% floor balances these sources.
                σ_W = max(precomp.σ_W_est * precomp.Wb_nodes[js[idx]], 0.15 * precomp.Wb_nodes[js[idx]])
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
                valid_ts = Int[], precomp = precomp,
                fallback = true)  # no inference possible = prior fallback
    end

    # Use data-estimated σ_obs if not explicitly provided
    σ_use = isnan(σ_obs) ? precomp.σ_obs_est : σ_obs

    # For data-sparse reaches, the likelihood has too little constraining
    # power for meaningful inference. Rather than producing unreliable
    # estimates, mark the reach as invalid and return fill values.
    n_valid = length(precomp.valid_ts)
    data_fraction = n_valid / reach.nt
    if data_fraction < 0.05
        @warn "infer: only $(n_valid)/$(reach.nt) valid timesteps ($(round(data_fraction*100, digits=1))%). Insufficient data for inference."
        nt = reach.nt
        return (Q_post  = fill(NaN, nt),
                Q_std    = fill(NaN, nt),
                n_post   = NaN, r_post = NaN, z0_post = NaN,
                θ_map    = fill(NaN, nt + 3),
                Σ        = zeros(nt + 3, nt + 3),
                result   = nothing,
                valid_ts = precomp.valid_ts,
                precomp  = precomp,
                fallback = true)
    end

    # Adaptive smoothness: reduce λ_smooth for data-sparse reaches so that
    # monthly priors can drive the seasonal Q cycle. When only a small
    # fraction of timesteps have observations, the smoothness penalty should
    # be relaxed to let the prior dominate the temporal structure.
    # Scale λ_smooth by data fraction: when data_fraction < 0.1, almost no
    # smoothness penalty (let monthly priors drive Q). When data_fraction > 0.5,
    # use the full λ_smooth (data has strong temporal structure).
    # This prevents the smoothness penalty from flattening Q on reaches with
    # very sparse observations.
    λ_use = if data_fraction < 0.1
        λ_smooth * data_fraction / 0.1  # linear ramp from 0 to λ_smooth
    else
        λ_smooth
    end
    println("  σ_obs = $(round(σ_use, digits=3)) m, λ_smooth = $(round(λ_use, digits=4)) (data_fraction=$(round(data_fraction, digits=3)))")

    # Initialize from prior medians
    θ0 = initialize_theta(p, precomp, months)

    # Objective closure
    obj(θ) = neg_log_joint(θ, precomp, p; σ_obs=σ_use, ν=ν, λ_smooth=λ_use, use_width=use_width)

    # In-place gradient via ForwardDiff
    function g!(G, θ)
        G .= ForwardDiff.gradient(obj, θ)
    end

    # L-BFGS optimization
    nt = reach.nt
    σ_obs_initial = σ_use   # save initial σ_obs to detect retries
    result = optimize(obj, g!, θ0, LBFGS(),
                      Optim.Options(iterations=iterations,
                                     g_tol=g_tol,
                                     allow_f_increases=true,
                                     show_trace=false))

    # Check for unphysical parameter values. If the optimizer produced
    # unphysical parameters (n outside [0.01, 0.1], r > 15, z0 far from hmin),
    # progressively increase σ_obs to allow more model–data mismatch.
    # This handles cases where σ_obs is too tight for small rivers where the
    # Manning + Dingman model has large structural error.
    #
    # NOTE: The previous logic only retried when !converged && !physical,
    # but convergence to unphysical values is equally problematic — it means
    # the prior was overwhelmed by the likelihood. We now also retry when
    # the optimizer converges to unphysical parameters.
    n_post_check  = exp(Optim.minimizer(result)[nt + 1])
    r_post_check  = exp(Optim.minimizer(result)[nt + 2])
    z0_post_check = Optim.minimizer(result)[nt + 3]

    # Physical bounds for static parameters.  These are slightly wider than
    # the prior bounds to allow the optimizer some room, but values well
    # outside indicate a degenerate solution.
    n_lo, n_hi = 0.005, 0.12
    r_hi = 20.0
    z0_lo, z0_hi = reach.hmin - 30, reach.hmin + 10

    physical = n_lo <= n_post_check <= n_hi &&
                r_post_check <= r_hi &&
                z0_lo <= z0_post_check <= z0_hi

    # Additional checks for inference quality:
    # 1. n well outside the prior support (n > 0.10 is physically unreasonable
    #    for any natural channel — the 0.05 truncation is regularisation, but
    #    values beyond 0.10 indicate a degenerate solution).
    # 2. Q values far outside the prior range (spike detection).
    n_outside_prior = n_post_check > 0.10

    # Quick Q sanity check: are estimated Q values wildly outside the prior?
    Q_check = exp.(Optim.minimizer(result)[1:nt])
    q_prior_medians = Float64[]
    for t in 1:nt
        ki = findfirst(==(t), precomp.valid_ts)
        mo = isnothing(ki) ? 1 : precomp.months_v[ki]
        push!(q_prior_medians, quantile(monthly_q_prior(p, mo), 0.5))
    end
    # Flag if any Q exceeds 10× the prior median (indicates unrealistic discharge)
    q_spike = any(Q_check .> 10.0 .* q_prior_medians)

    if !physical || n_outside_prior || q_spike
        @warn "infer: unphysical parameters (n=$(round(n_post_check,digits=4)), r=$(round(r_post_check,digits=2)), z0=$(round(z0_post_check,digits=2))), converged=$(Optim.converged(result)). Progressively increasing σ_obs."
        max_σ = 10.0 * σ_use
        while σ_use < max_σ
            σ_use *= 2.0
            obj_retry(θ) = neg_log_joint(θ, precomp, p; σ_obs=σ_use, ν=ν, λ_smooth=λ_use, use_width=use_width)
            function g_retry!(G, θ)
                G .= ForwardDiff.gradient(obj_retry, θ)
            end
            result = optimize(obj_retry, g_retry!, θ0, LBFGS(),
                              Optim.Options(iterations=iterations,
                                             g_tol=g_tol,
                                             allow_f_increases=true,
                                             show_trace=false))
            n_post_check  = exp(Optim.minimizer(result)[nt + 1])
            r_post_check  = exp(Optim.minimizer(result)[nt + 2])
            z0_post_check = Optim.minimizer(result)[nt + 3]
            physical = n_lo <= n_post_check <= n_hi &&
                        r_post_check <= r_hi &&
                        z0_lo <= z0_post_check <= z0_hi
            # Re-check Q sanity with updated parameters
            Q_retry = exp.(Optim.minimizer(result)[1:nt])
            q_spike_retry = any(Q_retry .> 10.0 .* q_prior_medians)
            n_outside_prior_retry = n_post_check > 0.10
            if physical && !n_outside_prior_retry && !q_spike_retry
                println("  Converged with σ_obs = $(round(σ_use, digits=3)) m")
                break
            end
        end
    end

    # Final check: are the parameters still unphysical or producing spikes after retries?
    Q_final = exp.(Optim.minimizer(result)[1:nt])
    q_spike_final = any(Q_final .> 10.0 .* q_prior_medians)
    final_physical = n_lo <= n_post_check <= n_hi &&
                      r_post_check <= r_hi &&
                      z0_lo <= z0_post_check <= z0_hi &&
                      n_post_check <= 0.10 &&
                      !q_spike_final

    # Additional fallback: when σ_obs was increased via retry AND the optimizer
    # pushes n beyond the prior truncation (0.05), the Manning+Dingman model
    # is compensating for structural error.  The WSE anomaly method is more
    # reliable in this case — it avoids overfitting the WSE data with an
    # excessively rough n that inflates discharge estimates.
    n_beyond_prior = n_post_check > 0.05
    σ_was_increased = σ_use > σ_obs_initial
    if !final_physical || (n_beyond_prior && σ_was_increased)
        if !final_physical
            @warn "infer: inference did not produce physical parameters after retries. Falling back to prior + WSE anomaly."
        else
            @warn "infer: n=$(round(n_post_check,digits=4)) exceeds prior truncation (0.05) after σ_obs retry — model compensating for structural error. Falling back to prior + WSE anomaly."
        end
        # --- Prior + WSE anomaly fallback ---
        # Use prior medians for n, r, z0 and compute Q from the WSE observations
        # using a simple depth-scaling anomaly relative to the prior.
        return _prior_anomaly_fallback(p, precomp, reach)
    end

    θ_map = Optim.minimizer(result)

    # Extract physical parameters
    Q_post  = exp.(θ_map[1:nt])
    n_post  = exp(θ_map[nt + 1])
    r_post  = exp(θ_map[nt + 2])
    z0_post = θ_map[nt + 3]

    # --- Profile refit: fix n, r and re-estimate Q, z₀ only ---
    # The joint MAP estimation can produce biased Q because the optimizer
    # compensates between Q and n (the Q–n degeneracy). By fixing n and r
    # at their MAP estimates and re-optimizing only Q and z₀, the WSE
    # observations directly constrain Q magnitude without the degeneracy.
    #
    # This works because n and r are primarily determined by the WSE profile
    # *shape* (curvature, backwater) rather than its absolute level. Once
    # shape parameters are fixed, Q and z₀ are the only degrees of freedom
    # affecting the WSE level, making Q identifiable from WSE alone.
    θ_profile = _profile_refit(θ_map, precomp, p; σ_obs=σ_use, ν=ν,
                                λ_smooth=λ_use, use_width=use_width)
    if θ_profile !== nothing
        Q_post  = exp.(θ_profile[1:nt])
        z0_post = θ_profile[nt + 1]  # z₀ may shift to accommodate new Q
        # n and r remain at their MAP values
        println("  Profile refit: Q range=$(round(minimum(Q_post),digits=1))..$(round(maximum(Q_post),digits=1)), z0=$(round(z0_post,digits=2))")
    end

    # Laplace uncertainty (use the final objective with the correct σ_obs)
    final_obj(θ) = neg_log_joint(θ, precomp, p; σ_obs=σ_use, ν=ν, λ_smooth=λ_use, use_width=use_width)
    Σ = laplace_uncertainty(final_obj, θ_map)
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
            precomp  = precomp,
            fallback = false)
end

"""
    _profile_refit(θ_map, precomp, priors; kwargs...) -> Vector or nothing

Profile refit: fix n and r at their MAP estimates and re-optimize Q and z₀ only.

After the joint MAP estimation (which estimates Q, n, r, z₀ simultaneously),
the n–Q degeneracy can cause Q to be biased toward the prior. By fixing
n and r (which are primarily determined by the WSE profile shape) and
re-optimizing only Q_{1:T} and z₀ (which determine the WSE level),
the WSE observations directly constrain Q magnitude.

The parameter vector for the profile refit is:
- θ[1:nt]      = log(Q_t)   — discharge at each timestep
- θ[nt + 1]  = z₀          — downstream bed elevation
- n and r are fixed at their MAP values from the joint estimation

# Returns
Optimized parameter vector of length nt+1, or `nothing` if the refit fails.
"""
function _profile_refit(θ_map, precomp::ManningPrecomp, priors::SWOTPriors;
                        σ_obs, ν, λ_smooth, use_width)
    nt = precomp.nt

    # Fix n and r at their MAP values
    n_fixed = exp(θ_map[nt + 1])
    r_fixed = exp(θ_map[nt + 2])

    # Yb depends on z₀, so we need to recompute it in the objective
    # Initialize from MAP: logQ from joint, z₀ from joint
    θ0_profile = vcat(θ_map[1:nt], θ_map[nt + 3])  # logQ_1..logQ_nt, z₀

    # Profile objective: fix n and r, optimize only logQ and z₀
    function profile_obj(θ_p)
        logQ_p = view(θ_p, 1:nt)
        z0_p   = θ_p[nt + 1]

        # Yb depends on z₀ — must recompute for each evaluation
        Yb_nodes = max.(precomp.hbf_z_nodes .- z0_p, 0.01)

        nll = zero(z0_p)  # correct Dual type

        # Pre-compute month index
        month_map = zeros(Int, nt)
        for (ki, t) in enumerate(precomp.valid_ts)
            month_map[t] = precomp.months_v[ki]
        end
        for t in 1:nt
            month_map[t] = month_map[t] == 0 ? 1 : month_map[t]
        end

        # --- 1. WSE likelihood (same as neg_log_joint but with n, r fixed) ---
        for k in 1:length(precomp.valid_ts)
            Q_t = exp(logQ_p[precomp.valid_ts[k]])
            WSE_pred = manning_wse_backwater(
                Q_t, n_fixed, r_fixed, z0_p, precomp.S0, precomp.H_bc[k],
                precomp.x_nodes, precomp.Wb_nodes, Yb_nodes, precomp.z_ref_nodes)

            js   = precomp.upstream_j[k]
            Hobs = precomp.upstream_H[k]
            Wobs = precomp.upstream_W[k]
            σ = isnan(σ_obs) ? precomp.σ_obs_est : σ_obs
            for idx in 1:length(js)
                if !isnan(Hobs[idx])
                    residual_H = Hobs[idx] - WSE_pred[js[idx]]
                    if isinf(ν)
                        nll += 0.5 * (residual_H / σ)^2
                    else
                        nll += ((ν + 1) / 2) * log(1 + (residual_H / σ)^2 / ν)
                end
            end
                if use_width && !isnan(Wobs[idx])
                    y_at_node = max((WSE_pred[js[idx]] - z0_p - precomp.z_ref_nodes[js[idx]]) * r_fixed / (r_fixed + 1), 0.01)
                    Yb_at_node = max(precomp.hbf_z_nodes[js[idx]] - z0_p, 0.01)
                    W_pred = precomp.Wb_nodes[js[idx]] * (y_at_node / Yb_at_node)^(1/r_fixed)
                    σ_W = max(precomp.σ_W_est * precomp.Wb_nodes[js[idx]], 0.15 * precomp.Wb_nodes[js[idx]])
                    residual_W = Wobs[idx] - W_pred
                    nll += 0.5 * (residual_W / σ_W)^2
                end
            end
        end

        # --- 2. Q priors (log-space, same as joint) ---
        for t in 1:nt
            q_prior = monthly_q_prior(priors, month_map[t])
            d = q_prior isa Truncated ? q_prior.untruncated : q_prior
            if d isa LogNormal
                nll -= logpdf(Normal(d.μ, d.σ), logQ_p[t])
            else
                nll -= safe_logpdf(q_prior, exp(logQ_p[t])) + logQ_p[t]
            end
        end

        # --- 3. z₀ prior ---
        nll -= safe_logpdf(priors.zp, z0_p)

        # --- 4. Temporal smoothness (same as joint) ---
        for t in 2:nt
            nll += λ_smooth * (logQ_p[t] - logQ_p[t - 1])^2
        end

        return nll
    end

    # Gradient via ForwardDiff
    function profile_g!(G, θ_p)
        G .= ForwardDiff.gradient(profile_obj, θ_p)
    end

    # L-BFGS optimization
    try
        result_p = optimize(profile_obj, profile_g!, θ0_profile, LBFGS(),
                           Optim.Options(iterations=500, g_tol=1e-6,
                                         allow_f_increases=true))
        if Optim.converged(result_p) || result_p.iterations > 10
            return Optim.minimizer(result_p)
        else
            return nothing
        end
    catch e
        @warn "_profile_refit: optimization failed" exception=(e, catch_backtrace())
        return nothing
    end
end

"""
    _prior_anomaly_fallback(p, precomp, reach) -> NamedTuple

Compute discharge using prior medians for channel parameters (n, r, z₀)
and a simple WSE-based scaling anomaly, for cases where MAP inference
fails to produce physical parameters.

The discharge for each timestep is estimated as:

    Q_t = Q_prior(t) × α_t

where `Q_prior(t)` is the monthly climatological prior median and `α_t`
is a scaling anomaly derived from the observed WSE change relative to what
the prior predicts.

The anomaly α_t is computed from the downstream boundary condition depth
and the prior-predicted uniform flow depth:

    α_t = (y_bc_t / y0_prior) ^ (5/3 + 1/r_prior)

where y_bc_t = (H_bc_t - z₀) × r/(r+1) is the observed downstream depth
(using the prior z₀) and y0_prior is the uniform flow depth from the prior
Q, n, and r. This exploits the Manning power-law Q ∝ y^(5/3+1/r).

Uncertainty is set to the prior standard deviation, reflecting that this
is a prior-informed estimate, not a data-constrained inference.
"""
function _prior_anomaly_fallback(p::SWOTPriors, precomp::ManningPrecomp,
                                   reach::SWOTReach)
    nt = precomp.nt

    # Use prior medians for static parameters
    n_prior = quantile(p.np, 0.5)
    r_prior = quantile(p.rp, 0.5)
    z0_prior = quantile(p.zp, 0.5)

    # Reach-averaged geometry
    Wb_mean = mean(precomp.Wb_nodes)
    Yb_prior = max.(precomp.hbf_z_nodes .- z0_prior, 0.01)
    Yb_mean = mean(Yb_prior)

    # Prior Q for each timestep (monthly)
    Q_prior_vec = fill(NaN, nt)
    Q_std_prior = fill(NaN, nt)
    for t in 1:nt
        ki = findfirst(==(t), precomp.valid_ts)
        mo = isnothing(ki) ? 1 : precomp.months_v[ki]
        q_prior = monthly_q_prior(p, mo)
        Q_prior_vec[t] = quantile(q_prior, 0.5)
        # Estimate prior standard deviation from quantiles
        # For LogNormal: std(Q) ≈ Q_median * (exp(σ²) - 1)^(1/2)
        # For truncated distributions, use interquartile range as a robust estimate
        q25 = quantile(q_prior, 0.25)
        q75 = quantile(q_prior, 0.75)
        # IQR-based std: std ≈ IQR / 1.35 (for Normal)
        Q_std_prior[t] = (q75 - q25) / 1.35
        if isnan(Q_std_prior[t]) || Q_std_prior[t] <= 0
            # Fall back to CV ≈ 1.0 for wide log-normal priors
            Q_std_prior[t] = Q_prior_vec[t]
        end
    end

    # Compute the uniform-flow depth for the prior Q
    # Manning depth: y = (Q*n / (C * √S))^(1/e), e = 5/3 + 1/r
    e_prior = 5.0/3.0 + 1.0/r_prior
    C_prior = Wb_mean * ((r_prior + 1) / (r_prior * Yb_mean))^(1/r_prior)

    Q_post = fill(NaN, nt)
    Q_std  = fill(NaN, nt)

    for (ki, t) in enumerate(precomp.valid_ts)
        Q_p = Q_prior_vec[t]

        # Uniform-flow depth for prior Q
        y0_prior = manning_depth(Q_p, n_prior, r_prior, Wb_mean, Yb_mean, precomp.S0)
        y0_prior = max(y0_prior, 0.01)

        # Observed downstream depth (boundary condition)
        y_bc = max((precomp.H_bc[ki] - z0_prior) * r_prior / (r_prior + 1), 0.01)

        # WSE anomaly: fractional change in depth relative to uniform flow.
        # Use a blended approach that stays close to the prior for small anomalies
        # and attenuates large anomalies to avoid extreme Q departures.
        #
        # The raw Manning power-law gives: Q/Q0 ≈ (y_bc/y0)^e
        # But this is highly sensitive to both model andprior errors.
        # Instead, use a log-linear anomaly:
        #   Q = Q_prior * exp(β * log(y_bc / y0))
        # with β < e to dampen the response and avoid runaway amplification.
        # β = 1.0 gives Q ∝ (y/y0), β = e gives the full nonlinear response.
        #
        # We use β = 1.5 as a conservative compromise between linearity (β=1)
        # and full Manning scaling (β≈2.5 for r≈3).
        depth_ratio = y_bc / y0_prior
        log_anomaly = 1.5 * log(max(depth_ratio, 0.1))  # β = 1.5

        # Blend: anomaly is attenuated toward zero (i.e., toward the prior)
        # The blend weight is σ = 0.5, meaning we only go 50% toward the anomaly.
        # This acknowledges that with poor inference, we should stay close to prior.
        α = exp(0.5 * log_anomaly)  # 50% blend toward prior

        # Hard limits: Q must stay within [0.2×, 5×] prior for safety
        α = clamp(α, 0.2, 5.0)

        Q_post[t] = Q_p * α

        # Uncertainty: prior uncertainty scaled by departure from prior
        # but at least as large as the prior std
        Q_std[t] = max(Q_std_prior[t], 0.3 * Q_post[t])
    end

    # For invalid timesteps, use prior median
    for t in 1:nt
        if isnan(Q_post[t])
            Q_post[t] = Q_prior_vec[t]
            Q_std[t] = Q_std_prior[t]
        end
    end

    println("  Fallback: using prior + WSE anomaly (n=$(round(n_prior,digits=4)), r=$(round(r_prior,digits=2)), z0=$(round(z0_prior,digits=2)))")
    println("  Q range: $(round(minimum(Q_post[.!isnan.(Q_post)]), digits=2)) .. $(round(maximum(Q_post[.!isnan.(Q_post)]), digits=2))")

    # Return with fallback flag
    return (Q_post  = Q_post,
            Q_std   = Q_std,
            n_post  = n_prior,
            r_post  = r_prior,
            z0_post = z0_prior,
            θ_map   = fill(NaN, nt + 3),  # no meaningful MAP vector
            Σ       = Diagonal(vcat(fill(NaN, nt), fill(NaN, 3))) |> Matrix,
            result  = nothing,
            valid_ts = precomp.valid_ts,
            precomp  = precomp,
            fallback = true)
end

"""
    initialize_theta(priors, precomp, months) -> Vector{Float64}

Construct initial parameter vector θ₀ from prior medians.

For Q, uses monthly climatological log-normal prior medians in log-Q space.
The `months` vector (length nt) maps each timestep to its calendar month,
allowing correct seasonal initialization even for unobserved timesteps.
"""
function initialize_theta(p::SWOTPriors, precomp::ManningPrecomp,
                           months::Vector{Int} = fill(1, precomp.nt))
    nt = precomp.nt

    # Static parameters from prior medians
    n_init = quantile(p.np, 0.5)
    r_init = quantile(p.rp, 0.5)
    z0_init = quantile(p.zp, 0.5)

    # Initialize Q from monthly prior medians using the correct calendar month
    # for each timestep. Previously, unobserved timesteps fell back to month 1,
    # which flattened the seasonal cycle for data-sparse reaches.
    logQ_init = fill(0.0, nt)
    for t in 1:nt
        mo = months[t]
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
function compute_A0(reach::SWOTReach, res::NamedTuple)
    r_mean = res.r_post
    z0_mean = res.z0_post

    # downstream node geometry
    x_ds = reach.x[1]
    Wb   = reach.wbf(x_ds)                      # bankfull width [m]
    hbf  = reach.hbf(x_ds)                      # bankfull WSE [m]
    zx   = z0_mean + reach.z(x_ds)              # bed elevation at x_ds
                                                  # = z0_mean since z(x[1])=0
    Yb   = max(hbf - zx, 0.01)                  # bankfull mean depth [m]

    # mean flow depth at minimum observed WSE
    y0   = max((reach.hmin - zx) * r_mean / (r_mean + 1), 0.01)

    # Dingman area at y0
    Ym   = (r_mean + 1) / r_mean * y0
    A0   = Wb * (Ym / Yb)^(1 / r_mean) * y0

    return A0
end
