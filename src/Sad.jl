module Sad

abstract type CrossSection end

include("crosssections.jl")
include("gvf.jl")
include("priors.jl")
include("kalman.jl")

"""
    assimilate(H::Vector{Float64})

Assimilate SWOT observations for river reach.

# Arguments

- `H`: water surface elevation
- `W`: water surface width
- `x`: downstream distance for each cross section
- `wbf`: bankfull width
- `hbf`: bankfull depth
- `Qₚ`: prior probability distribution for discharge
- `nₚ`: prior probability distribution for roughness coefficient
- `rₚ`: prior probability distribution for channel shape parameter
- `zₚ`: prior distribution for downstream bed elevation
- `nens`: ensemble size
- `ϵₒ`: observation error standard deviation

"""
function assimilate(H, W, x, wbf, hbf, S, Qp, np, rp, zp, nens, ϵ₀::Float64=0.01)
    seed!(1)
    min_ensemble_size = 5
    nt = size(H, 2)
    Qa = zeros(1, nt)
    na = zeros(1, nt)
    # FIXME: estimate bed slope
    S0 = S
    Qe, ne, re, ze = prior_ensemble(x, Qp, np, rp, zp, nens)
    for t in 1:nt
        he = gvf_ensemble!(H[1, t], S0, x, hbf, wbf, Qe, ne, re, ze)
        # find valid flow depths simulated
        i = findall(he[1, :] .> 0)
        if length(i) < min_ensemble_size
            Sf = [S[j, t] > 0 ? S[j, t] : minimum(S[:, t][S[:, t] .> 0]) for j=1:length(x)]
            he = [(Qe[e] .* ne[e]) ./ (W[j, t] .* Sf[j].^0.5).^(3/5) for j=1:length(x), e=1:nens]
        end
        i = findall(he[1, :] .> 0)
        h = ze .+ he .* ((re .+ 1) ./ re)'
        X = zeros(2, length(i))
        X[1, :] = Qe[i]
        X[2, :] = ne[i]
        XA = h[:, i]
        d = H[:, 1]
        E = rand(Normal(ϵ₀, ϵ₀/1000), length(d), length(i)) .* rand([-1, 1], length(d), length(i))
        A = letkf(X, d, XA, E, diagR=true)
        Qa[t] = mean(A[1, :])
        na[t] = mean(A[2, :])
    end
    Qa, na
end

"""
    bathymetry!(ze, Se, Qe, ne, re, x, hbf, wbf, H, ϵ₀::Float64=0.01)

Estimate channel bed slope and thalweg (elevation) by assimilating SWOT water surface elevation observations.

# Arguments
- `ze`: ensemble of bed elevations
- `Se`: ensemble of bed slopes
- `Qe`: ensemble of discharge
- `ne`: ensemble of roughness coefficient
- `re`: ensemble of channel shape parameters
- `x`: channel chainage
- `hbf`: bankfull water surface elevation
- `wbf`: bankfull width
- `H`: time-averaged water surface elevation observed
- `ϵ₀`: observation error standard deviation

"""
function bathymetry!(ze::Matrix{Float64}, Se::Matrix{Float64}, Qe::Vector{Float64}, ne::Vector{Float64}, re::Vector{Float64}, x::Vector{Float64}, hbf::Vector{Float64}, wbf::Vector{Float64}, H::Vector{Float64}, ϵ₀::Float64=0.01)
    seed!(1)
    min_ensemble_size = 5
    he = gvf_ensemble!(H[1], Se, x, hbf, wbf, Qe, ne, re, ze)
    i = findall(he[1, :] .> 0)
    if length(i) > min_ensemble_size
        h = ze .+ he .* ((re .+ 1) ./ re)'
        X = ze[:, i]
        XA = h[:, i]
        d = H
        E = rand(Normal(ϵ₀, ϵ₀/1000), length(d), length(i)) .* rand([-1, 1], length(d), length(i))
        A = letkf(X, d, XA, E, diagR=true)
        # FIXME ensure that estimated bathymetry is plausible
        ze = A
        Se = diff(ze, dims=1) ./ diff(x)
    end
    nothing
end

"""
    estimate(x, H, W, Qp, np, rp, zp, nens)

Estimate discharge and flow parameters from SWOT observations.

# Arguments
- `x`: channel chainage
- `H`: time series of water surface elevation profiles
- `W`:time series of width profiles
- `Qp`: discharge prior distribution
- `np`: roughness coefficient prior distribution
- `rp`: channel shape parameter prior distribution
- `zp`: downstream bed elevation prior distribution
- `nens`: ensemble size

"""
function estimate(x::Vector{Float64}, H::Matrix{Float64}, W::Matrix{Float64}, Qp::Distribution, np::Distribution, rp::Distribution, zp::Distribution, nens::Int)
    wbf = maximum(W, dims=2)[:, 1]
    hbf = maximum(H, dims=2)[:, 1]
    S0 = mean(diff(H, dims=1) ./ diff(x), dims=2)[:, 1]
    S0 = [S0[1]; S0]
    Qp, np, rp, zp = rejection_sampling(Qp, np, rp, zp, x, H, S0, mean(H[1, :]), wbf, hbf, nens);
    Qe, ne, re, ze = prior_ensemble(x, Qp, np, rp, zp, nens);
    S = diff(H, dims=1) ./ diff(x);
    Se = repeat(S', outer=((nens ÷ size(H, 2)) + 1))'[:, 1:nens];
    Se = [Se[1, :]'; Se]
    bathymetry!(ze, Se, Qe, ne, re, x, hbf, wbf, mean(H, dims=2)[:, 1])
    zp = Truncated(Normal(mean(ze[1, :]), std(ze[1, :])), -Inf, minimum(H[1, :]))
    Qa, na = assimilate(H, W, x, wbf, hbf, Se, Qp, np, rp, zp, nens)
    Qa, na
end

export
    width,
    depth,
    gvf,
    assimilate


end
