using Distributions
using Random: seed!
using KernelDensity: kde
using BlackBoxOptim.Utils: latin_hypercube_sampling

"""
    get_samples(p, samples)

Get quantiles from `p` according to weights provided by
Latin Hypercube Sampling `samples`.

"""
function get_samples(p, samples)
    out = try
        quantile.(p, samples)
    catch
        [quantile(p, s) for s in samples]
    end
    out
end

"""
    lhs_ensemble(nens, args...)

Generate an ensemble of size `nens` using Latin Hypercube Sampling of the
list of distributions or collections provided as `args`.

"""
function lhs_ensemble(nens, args...)
    seed!(1)
    nvars = length(args)
    usamples = latin_hypercube_sampling([0. for _ in 1:nvars], [1. for _ in 1:nvars], nens)
    [get_samples(args[i], usamples[i, :]) for i in 1:nvars]
end

"""
    prior_ensemble(x::Vector{Float64}, Qp::Distribution, np::Distribution, rp::Distribution, zp::Distribution, nens::Int)

Generate a prior ensemble of discharge, roughness coefficient, channel shape parameter,
and downstream bed elevation from provided distributions or sample data.

# Arguments
- `x`: distances between cross sections
- `Qp`: discharge prior distribution
- `np`: roughness coefficient distribution
- `rp`: channel shape parameter distribution
- `zp`: downstream bathymetry distribution
- `nens`: ensemble size

"""
function prior_ensemble(x::Vector{Float64}, Qp::Distribution, np::Distribution,
                        rp::Distribution, zp::Distribution, nens::Int)
    ze = zeros(length(x), nens)
    Qe, ne, re, ze[1, :] = lhs_ensemble(nens, Qp, np, rp, zp)
    Qe, ne, re, ze
end

"""
    gvf_ensemble!(H::Vector{Float64}, S, x::Vector{Float64}, hbf::Vector{Float64}, wbf::Vector{Float64}, Qe::Vector{Float64}, ne::Vector{Float64}, re::Vector{Float64}, ze)

Generate an ensemble of water height profiles from Gradually-Varied-Flow simulations
and associated profiles of bed elevation.

- `H`: water surface elevation
- `S`: bed slope
- `x`: channel chainage
- `hbf`: bankfull water surface elevation
- `wbf`: bankfull width
- `Qe`: ensemble discharge
- `ne`: ensemble roughness coefficient
- `re`: ensemble channel shape parameter
- `ze`: ensemble bed elevation profiles

"""
function gvf_ensemble!(hbc::Float64, S, x::Vector{Float64}, hbf::Vector{Float64}, wbf::Vector{Float64},
                       Qe::Vector{Float64}, ne::Vector{Float64}, re::Vector{Float64}, ze)
    nens = length(Qe)
    if ndims(S) > 1
        Se = S
    else
        Se = repeat(S', outer=nens)'
    end
    for j in 2:length(x)
        ze[j, :] = ze[j-1, :] .+ Se[j, :] .* (x[j] - x[j-1]);
    end
    ybf = hbf .- ze
    he = zeros(length(x), nens)
    for i in 1:nens
        he[:, i] = gvf(Qe[i], (hbc-ze[1, i])*re[i]/(re[i] + 1), Se[:, i], ne[i], x, wbf, ybf[:, i], [re[i] for _ in 1:length(x)])
    end
    he
end

"""
    rejection_sampling(Qp, np, rp, zp, x, H, S0, hbc, wbf, hbf, nens)

Use rejection sampling to select a subset of the prior ensemble.

# Arguments
- `Qp`: prior discharge distribution
- `np`: prior roughness coefficient distribution
- `rp`: prior channel shape parameter distribution
- `zp`: prior distribution for downstream bed elevation
- `x`: channel chainage
- `H`: time series of observed water surface elevation profiles
- `S0`: prior estimate of channel bed slope
- `hbc`: mean downstream flow depth used as boundary condition
- `wbf`: bankfull width
- `hbf`: bankfull water surface elevation
- `nens`: ensemble size

"""
function rejection_sampling(Qp::Distribution, np::Distribution, rp::Distribution, zp::Distribution,
                            x::Vector{Float64}, H::Matrix{Float64}, S0::Vector{Float64}, hbc::Float64,
                             wbf::Vector{Float64}, hbf::Vector{Float64}, nens::Int)
    Qe, ne, re, ze = prior_ensemble(x, Qp, np, rp, zp, nens)
    he = gvf_ensemble!(hbc, S0, x, hbf, wbf, Qe, ne, re, ze)
    h = ze .+ he .* ((re .+ 1) ./ re)'
    i = findall(he[1, :] .> 0)
    obs = mean(H, dims=1)[1, :]
    Fobs = kde(obs)
    mod = mean(h[:, i], dims=1)[1, :]
    Fmod = kde(mod)
    L = 1.0
    accepted = [rand(Uniform(0, L)) * pdf(Fmod, s) <= pdf(Fobs, s) for s in mod]
    Qpₘ = LogNormal(log(mean(Qe[i[accepted]])) - (std(Qe[i[accepted]]) / mean(Qe[i[accepted]]))^2/2,
                    std(Qe[i[accepted]]) / mean(Qe[i[accepted]]))
    rpₘ = Truncated(Normal(mean(re[i[accepted]]), std(re[i[accepted]])), 0.5, 20.0)
    zpₘ = Truncated(Normal(mean(ze[1, i[accepted]]), std(ze[1, i[accepted]])), -Inf, minimum(H[1, :]))
    Qpₘ, np, rpₘ, zpₘ
end

"""
    priors(ncfile)

Derive prior distributions for discharge, roughness coefficient, channel shape parameter, and bed elevation by reading the [SWORD](http://gaia.geosci.unc.edu/SWORD/) a-priori database.

# Arguments
- `ncfile`: NetCDF file of SWORD database

"""
function priors(ncfile::String)
    nothing
end

"""
    priors(qwbm, hmin)

Derive distributions for discharge, roughness coefficient, channel shape parameter, and bed elevation using uninformative priors.

# Arguments
- `qwbm`: mean discharge
- `hmin`: minimum downstream water surface elevation

"""
function priors(qwbm::Float64, hmin::Float64)
    Qp = Truncated(LogNormal(log(qwbm)-2.0^2/2, 2.0), 0.1*qwbm, 10*qwbm)
    np = Uniform(0.01, 0.07)
    rp = Truncated(Normal(2.5, 0.5), 0.5, 20.0)
    zp = Uniform(hmin-20, hmin)
    Qp, np, rp, zp
end
