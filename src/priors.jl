using Distributions
using NCDatasets
using Dates

"""
River planform classification used to set uninformative prior bounds on
the Dingman shape exponent `r`.
"""
@enum River braided=1 sinuous=2 more_sinuous=3 straight=4

"""
    SWOTPriors

Prior distributions for all inferred parameters in the SAD algorithm.

# Fields
- `Qp`: discharge priors — either a single `Distribution` (time-invariant)
        or a `Vector{Distribution}` of length 12 (one per calendar month)
        when monthly ML priors are available from the SoS database.
- `np`: Manning roughness coefficient prior
- `rp`: Dingman shape exponent prior
- `zp`: downstream bed elevation prior
"""
struct SWOTPriors
    Qp :: Union{Distribution, Vector{<:Distribution}}
    np :: Distribution
    rp :: Distribution
    zp :: Distribution
end

"""
    monthly_q_prior(priors, month) -> Distribution

Return the discharge prior for a given calendar month (1-12).
If `priors.Qp` is a single distribution (no monthly data), returns it
directly. If it is a 12-element vector, returns the distribution for
the given month.
"""
function monthly_q_prior(priors::SWOTPriors, month::Int)
    if priors.Qp isa Vector
        return priors.Qp[month]
    else
        return priors.Qp
    end
end

"""
    priors(sosfile, hmin, reachid)

Derive prior distributions from the SoS (SWORD of Science) database.

Uses monthly ML discharge priors if available (`monthly_q` group),
otherwise falls back to the time-invariant model mean.

# Arguments
- `sosfile`: path to SoS NetCDF file
- `hmin`:    minimum observed downstream WSE (upper bound for z0 prior)
- `reachid`: SWORD reach ID

Returns a `SWOTPriors` or `missing` if discharge prior cannot be constructed.
"""
function priors(sosfile::String, hmin::Float64, reachid::Int)
    Dataset(sosfile) do f
        gr   = f.group["reaches"]
        idxs = findall(gr["reach_id"][:] .== reachid)
        isempty(idxs) && return missing
        i = idxs[1]
        # roughness
        g = f.group["gbpriors"].group["reach"]
        n_l = exp(g["lowerbound_logn"][i])
        n_u = exp(g["upperbound_logn"][i])
        np = try Uniform(n_l, n_u) catch; Uniform(0.01, 0.07) end
        # channel shape parameter
        r_m = exp(g["logr_hat"][i])
        r_s = exp(g["logr_sd"][i])
        r_l = exp(g["lowerbound_logr"][i])
        r_u = exp(g["upperbound_logr"][i])
        rp  = try truncated(Normal(r_m, r_s), r_l, r_u) catch; Uniform(0.5, 10.0) end
        # discharge
        m   = NCDatasets.group(f, "model")
        q_m = m["mean_q"][i]
        q_u = m["max_q"][i]
        q_l = m["min_q"][i]
        if ismissing(q_m)
            return missing
        end
        Qp = _build_q_prior(f, i, q_m, q_l, q_u)
        # bed elevation — depth estimate scales with discharge magnitude
        q_m_f     = Float64(q_m)
        depth_est = q_m_f > 500.0 ? 7.0 : (q_m_f > 100.0 ? 5.0 : 3.0)
        z0_est    = hmin - depth_est
        zp = Uniform(z0_est - 3.0, z0_est + 3.0)
        SWOTPriors(Qp, np, rp, zp)
    end
end

"""
    _build_q_prior(f, i, q_m, q_l, q_u) -> Distribution or Vector{Distribution}

Build discharge prior(s) from SoS dataset. Uses monthly ML priors if the
`monthly_q` group exists and has valid data for reach index `i`, otherwise
falls back to a single time-invariant truncated LogNormal.
"""
function _build_q_prior(f::NCDataset, i::Int,
                        q_m::Real, q_l::Real, q_u::Real)
    g = f.group["model"]
    if haskey(g, "monthly_q") && size(g["monthly_q"], 1) == 12
        monthly_means = g["monthly_q"][:, i]
        if !any(ismissing, monthly_means) && all(monthly_means .> 0)
            return _monthly_distributions(monthly_means, q_l, q_u, q_m)
        end
    end
    # fall back to time-invariant prior
    return _single_q_prior(q_m, q_l, q_u)
end

"""
    _monthly_distributions(monthly_means, q_l, q_u, q_m_annual)

Build a 12-element Vector{Distribution} of truncated LogNormals,
one per calendar month. Each distribution is centred on the monthly
ML mean with σ=2 in log-space (wide enough to cover uncertainty in
the ML estimate), truncated to [q_l, q_u * 20].
"""
function _monthly_distributions(monthly_means::Vector,
                                 q_l::Real, q_u::Real,
                                 q_m_annual::Real)
    map(1:12) do mo
        qm_mo = Float64(monthly_means[mo])
        logmu = log(qm_mo) - 1.0^2 / 2
        logmu = isinf(logmu) ? log(q_m_annual) : logmu
        q_hi  = max(q_u, 20 * qm_mo)
        q_lo  = max(q_l, 0.01)
        try
            truncated(LogNormal(logmu, 2.0), q_lo, q_hi)
        catch
            truncated(LogNormal(logmu, 2.0), 0.1 * qm_mo, 20 * qm_mo)
        end
    end
end

"""
    _single_q_prior(q_m, q_l, q_u) -> Distribution

Build a single time-invariant truncated LogNormal discharge prior.
"""
function _single_q_prior(q_m::Real, q_l::Real, q_u::Real)
    qm = log(q_m) - 2.0^2 / 2
    qm = isinf(qm) ? (q_u + q_l) / 2.0 : qm
    try
        truncated(LogNormal(qm, 2.0), q_l, q_u)
    catch
        truncated(LogNormal(qm, 2.0), 0.1 * q_m, 20 * q_m)
    end
end

"""
    priors(qwbm, hmin, class; reach)

Construct uninformative priors when SoS data are unavailable.

# Arguments
- `qwbm`:  mean discharge estimate [m³/s]
- `hmin`:  minimum observed downstream WSE [m]
- `class`: river planform class (`River` enum: braided=1, sinuous=2,
           more_sinuous=3, straight=4) — controls r prior bounds

# Keyword arguments
- `reach`: optional SWOTReach — if provided, z0 is estimated from the
           observed WSE profile and bed slope rather than hmin alone
"""
function priors(qwbm::Float64, hmin::Float64, class::River;
                reach=nothing)
    rbnds = [(0.5, 1.0), (1.0, 5.0), (5.0, 10.0), (10.0, 20.0)]
    q_cv = 1.0
    Qp = truncated(LogNormal(log(qwbm) - q_cv^2 / 2, q_cv), 0.1 * qwbm, 10 * qwbm)
    n_lo = qwbm > 500.0 ? 0.025 : (qwbm > 100.0 ? 0.020 : 0.015)
    np = Uniform(n_lo, 0.07)
    rp = Uniform(rbnds[Int(class)]...)
    zp = z0_prior(qwbm, hmin, reach)
    SWOTPriors(Qp, np, rp, zp)
end

"""
    z0_prior(qwbm, hmin, reach) -> Uniform

Estimate the downstream bed elevation prior. If reach is provided, uses
the minimum observed WSE minus an estimated depth from Manning scaling.
Otherwise falls back to a fixed depth estimate based on qwbm.
"""
function z0_prior(qwbm::Float64, hmin::Float64, reach)
    if !isnothing(reach)
        S_med = reach.S0(reach.x[end])
        W_med = reach.wbf(reach.x[end])
        n_est = qwbm > 500.0 ? 0.035 : 0.030
        S_use = max(S_med, 1e-5)
        depth_est = (n_est * qwbm / (W_med * sqrt(S_use)))^0.6
        depth_est = clamp(depth_est, 2.0, 20.0)
    else
        depth_est = qwbm > 500.0 ? 7.0 : (qwbm > 100.0 ? 5.0 : 3.0)
    end
    z0_est = hmin - depth_est
    return Uniform(z0_est - 3.0, z0_est + 3.0)
end
