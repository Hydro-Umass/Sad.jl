using Distributions
using NCDatasets

"""
River planform classification used to set uninformative prior bounds on
the Dingman shape exponent `r`.
"""
@enum River braided=1 sinuous=2 more_sinuous=3 straight=4

"""
    SWOTPriors

Prior distributions for all inferred parameters in the SAD algorithm.

# Fields
- `Qp`: discharge prior
- `np`: Manning roughness coefficient prior
- `rp`: Dingman shape exponent prior
- `zp`: downstream bed elevation prior
- `ap`: slope correction factor prior (centered on 1)
"""
struct SWOTPriors
    Qp :: Distribution
    np :: Distribution
    rp :: Distribution
    zp :: Distribution
    ap :: Distribution
end

"""
    priors(sosfile, hmin, reachid)

Derive prior distributions from the SoS (SWORD of Science) database.

# Arguments
- `sosfile`: path to SoS NetCDF file
- `hmin`:    minimum observed downstream WSE (upper bound for z0 prior)
- `reachid`: SWORD reach ID

Returns a `SWOTPriors` or `missing` if discharge prior cannot be constructed.
"""
function priors(sosfile::String, hmin::Float64, reachid::Int)
    NCDataset(sosfile) do f
        gr = NCDatasets.group(f, "reaches")
        i  = findall(gr["reach_id"][:] .== reachid)[1]
        # roughness
        g   = NCDatasets.group(NCDatasets.group(f, "gbpriors"), "reach")
        n_l = exp(g["lowerbound_logn"][i])
        n_u = exp(g["upperbound_logn"][i])
        np  = try Uniform(n_l, n_u) catch; Uniform(0.01, 0.10) end
        # channel shape parameter
        r_m = exp(g["logr_hat"][i])
        r_s = exp(g["logr_sd"][i])
        r_l = exp(g["lowerbound_logr"][i])
        r_u = exp(g["upperbound_logr"][i])
        rp  = try truncated(Normal(r_m, r_s), r_l, r_u) catch; Uniform(1.0, 10.0) end
        # discharge
        m   = NCDatasets.group(f, "model")
        q_m = m["mean_q"][i]
        q_u = m["max_q"][i]
        q_l = m["min_q"][i]
        if ismissing(q_m)
            return missing
        end
        qm = log(q_m) - 2.0^2 / 2
        qm = isinf(qm) ? (q_u + q_l) / 2.0 : qm
        # Upper bound widened to 20x mean Q to capture extreme flood events
        # beyond the 10x bound which excludes truth at high-flow reaches
        Qp = try
            truncated(LogNormal(qm, 2.0), q_l, q_u)
        catch
            truncated(LogNormal(qm, 2.0), 0.1 * q_m, 20 * q_m)
        end
        # bed elevation
        z0_est = hmin - 5.0
        zp = Uniform(z0_est - 3.0, z0_est + 3.0)
        # slope correction
        ap = LogNormal(0.0, 0.2)
        SWOTPriors(Qp, np, rp, zp, ap)
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
    # Upper bound widened to 20x mean Q to capture extreme flood events
    # (10x was excluding truth at high-flow reaches)
    Qp = truncated(LogNormal(log(qwbm) - 2.0^2 / 2, 2.0), 0.1 * qwbm, 20 * qwbm)
    n_lo = qwbm > 500.0 ? 0.025 : (qwbm > 100.0 ? 0.020 : 0.015)
    np = Uniform(n_lo, 0.07)
    rp = Uniform(rbnds[Int(class)]...)
    zp = z0_prior(qwbm, hmin, reach)
    ap = LogNormal(0.0, 0.2)
    SWOTPriors(Qp, np, rp, zp, ap)
end

"""
    z0_prior(qwbm, hmin, reach) -> Uniform

Estimate the downstream bed elevation prior. If reach is provided, uses
the minimum observed WSE minus an estimated depth from Manning scaling.
Otherwise falls back to a fixed depth estimate based on qwbm.
"""
function z0_prior(qwbm::Float64, hmin::Float64, reach)
    if !isnothing(reach)
        S_med = reach.S0(reach.x[end] / 2)
        W_med = reach.wbf(reach.x[end] / 2)
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
