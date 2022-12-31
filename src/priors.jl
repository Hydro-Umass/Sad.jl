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
    prior_ensemble(x, Qp, np, rp, zp, nens)

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
    rejection_sampling(Qp, np, rp, zp, x, H, S0, hbc, wbf, hbf, nens, nsamples)

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
- `nsamples`: number of samples

"""
function rejection_sampling(Qp::Distribution, np::Distribution, rp::Distribution, zp::Distribution, x::Vector{Float64}, H::Matrix{FloatM}, S0::Vector{Float64}, wbf::Vector{Float64}, hbf::Vector{Float64}, nens::Int, nsamples::Int)
    Qe, ne, re, ze = prior_ensemble(x, Qp, np, rp, zp, nsamples)
    he = gvf_ensemble!(mean(skipmissing(H[1, :])), S0, x, hbf, wbf, Qe, ne, re, ze)
    h = ze .+ he .* ((re .+ 1) ./ re)'
    i = findall(he[1, :] .> 0)
    # obs = mean(H, dims=1)[1, :]
    obs = mean.(skipmissing.(eachcol(H)))
    Fobs = kde(obs[.!isnan.(obs)])
    mod = mean(h[:, i], dims=1)[1, :]
    Fmod = kde(mod)
    L = 1.0
    accepted = [rand(Uniform(0, L)) * pdf(Fmod, s) <= pdf(Fobs, s) for s in mod]
    Qm = mean(Qe[i[accepted]])
    Qcv = std(Qe[i[accepted]]) / mean(Qe[i[accepted]])
    Qpₘ = truncated(LogNormal(log(Qm) - Qcv^2/2, Qcv), minimum(Qp), maximum(Qp))
    zpₘ = truncated(Normal(mean(ze[1, i[accepted]]), std(ze[1, i[accepted]])), -Inf, minimum(skipmissing(H[1, :])))
    re = lhs_ensemble(nens, rp)[1]
    Qe, ne, _, ze = prior_ensemble(x, Qpₘ, np, rp, zpₘ, size(H, 2))
    hm = zeros(nens)
    for e=1:nens
        he = zeros(length(x), size(H, 2))
        for t=1:size(H, 2)
            hbc = ismissing(H[1, t]) ? mean(skipmissing(H[1, :])) : H[1, t]
            he[:, t] = gvf_ensemble!(hbc, S0, x, hbf, wbf, [Qe[t]], [ne[t]], [re[e]], ze)
            h = ze .+ he .* ((re[e] .+ 1) ./ re[e])'
            i = findall(he[1, :] .> 0)
            hm[e] = mean(h[:, i])
        end
    end
    reₘ = re[argmin(abs.(mean(mean.(skipmissing.(eachrow(H)))) .- hm))]
    # assume a CV of 0.1 for the channel shape parameter
    rpₘ = truncated(Normal(reₘ, 0.1*reₘ), 0.5, 20.0)
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
function priors(qwbm::Float64, hmin::Float64, class::River)
    rbnds = [(0.5, 1), (1, 5), (5 ,10), (10, 20)]
    Qp = truncated(LogNormal(log(qwbm)-2.0^2/2, 2.0), 0.1*qwbm, 10*qwbm)
    np = Uniform(0.01, 0.07)
    # rp = truncated(Normal(2.5, 0.5), 0.5, 20.0)
    rp = Uniform(rbnds[Int(class)]...)
    zp = Uniform(hmin-20, hmin)
    Qp, np, rp, zp
end
