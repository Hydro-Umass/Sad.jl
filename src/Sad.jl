module Sad

abstract type CrossSection end

@enum River begin
    braided=1
    sinuous=2
    more_sinuous=3
    straight=4
end

include("crosssections.jl")
include("gvf.jl")
include("priors.jl")
include("kalman.jl")

"""
    assimilate(H, W, x, wbf, hbf, S0, Qp, np, rp, zp, nens, constrain)

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
- `constrain`: switch for applying AHG constraint

"""
function assimilate(H, W, x, wbf, hbf, S0::Vector{Float64}, Qp::Distribution, np::Distribution, rp::Distribution, zp::Distribution, nens::Int, constrain::Bool=true)
    seed!(1)
    min_ensemble_size = 5
    nt = size(H, 2)
    Qa = zeros(nt)
    Qu = zeros(nt)
    na = zeros(nt)
    S = diff(H, dims=1) ./ diff(x)
    S = [S[1, :]'; S]
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
        d = H[:, t]
        # adaptive estimation of observation error covariance
        E = (((d .- mean(XA,dims=2)).^2 .- std(XA, dims=2).^2))
        A = letkf(X, d, XA, E, diagR=true)
        # if estimated discharge is negative, use a log transform in the assimilation
        qm = mean(A[1, :])
        if qm < 0
            X[1, :] = log.(Qe[i])
            A = letkf(X, d, XA, E, diagR=true)
            A[1, :] .= exp.(A[1, :])
        end
        Aq = A[1, :]
        if constrain
            ahg_constrain!(Aq, H[1, t], hbf, wbf, S0, x, ne[i], re[i], ze[:, i])
        end
        Qa[t] = mean(Aq)
        Qu[t] = std(Aq)
        na[t] = mean(A[2, :])
    end
    Qa, Qu, na
end

"""
    bathymetry!(ze, Se, Qe, ne, re, x, hbf, wbf, H)

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

"""
function bathymetry!(ze::Matrix{Float64}, Se::Matrix{Float64}, Qe::Vector{Float64}, ne::Vector{Float64}, re::Vector{Float64}, x::Vector{Float64}, hbf::Vector{Float64}, wbf::Vector{Float64}, H::Vector{Float64})
    seed!(1)
    min_ensemble_size = 5
    he = gvf_ensemble!(H[1], Se, x, hbf, wbf, Qe, ne, re, ze)
    i = findall(he[1, :] .> 0)
    if length(i) > min_ensemble_size
        h = ze .+ he .* ((re .+ 1) ./ re)'
        X = ze[:, i]
        XA = h[:, i]
        d = H
        E = (((d .- mean(XA,dims=2)).^2 .- std(XA, dims=2).^2))
        A = letkf(X, d, XA, E, diagR=true)
        # FIXME ensure that estimated bathymetry is plausible
        ze = A
        Se = diff(ze, dims=1) ./ diff(x)
    end
    nothing
end

"""
    ahg_constrain!(Qa, hbc, hbf, wbf, S0, x, ne, re, ze, tol=1e-3)

Constrain assimilated discharge with At-a-Station hydraulic geometry (AHG) relationships.

# Arguments
- `Qa`: assimilated ensemble discharge
- `hbc`: downstream boundary water surface elevation
- `hbf`: bankfull water surface elevation
- `wbf`: bankfull width
- `S0`: bed slope
- `x`: channel chainage
- `ne`: ensemble roughness coefficient
- `re`: ensemble channel shape parameter
- `ze`: ensemble bed elevation
- `tol`: numeric tolerance for closing the AHG relationships

"""
function ahg_constrain!(Qa::Vector{Float64}, hbc::Float64, hbf::Vector{Float64}, wbf::Vector{Float64}, S0::Vector{Float64}, x::Vector{Float64}, ne::Vector{Float64}, re::Vector{Float64}, ze::Matrix{Float64}, tol::Float64=1e-3)
    j = findall(Qa .> 0)
    ha = gvf_ensemble!(hbc, S0, x, hbf, wbf, Qa[j], ne[j], re[j], ze[:, j])
    i = findall(ha[1, :] .> 0)
    ij = j[i] # index that combines positive discharge and flow depth
    w = zeros(length(x), length(i))
    for e=1:length(i)
        for k=1:length(x)
            w[k, e] = width(Dingman(wbf[k], hbf[k], ha[k, i[e]]*(re[ij[e]]+1)/re[ij[e]], re[ij[e]], S0[k], ne[ij[e]]))
        end
    end
    qₗ = log.(Qa[ij])
    wₗ = log.(w[1, :])
    yₗ = log.(ha[1, i] .* (re[ij] .+ 1) ./ re[ij])
    bw = (sum(qₗ .* wₗ) .- (1/length(qₗ)) .* sum(qₗ) .* sum(wₗ)) / (sum(qₗ.^2) .- (1/length(qₗ)) * sum(qₗ)^2)
    aw = mean(wₗ) - bw * mean(qₗ)
    by = (sum(qₗ .* yₗ) .- (1/length(qₗ)) .* sum(qₗ) .* sum(yₗ)) / (sum(qₗ.^2) .- (1/length(qₗ)) * sum(qₗ)^2)
    ay = mean(yₗ) - by * mean(qₗ)
    X = zeros(1, length(ij))
    X[1, :] = qₗ
    XA = zeros(2, length(ij))
    XA[1, :] = sqrt.((wₗ .- (aw .+ bw .* qₗ)).^2)
    XA[2, :] = sqrt.((yₗ .- (ay .+ by .* qₗ)).^2)
    d = [0.0; 0.0]
    E = diagm([tol, tol])
    A = letkf(X, d, XA, E, diagR=true)
    Qa .= ((Qa .- mean(Qa[j])) ./ std(Qa)) .* std(exp.(A)) .+ mean(exp.(A))
end

"""
    flow_parameters(Qa, na, x, H, W, S, S0, hbf, wbf, r, z)

Estimate flow parameters (roughness coefficient and baseflow cross-sectional area) from assimilated discharge.

# Arguments
- `Qa`: time series of estimated discharge
- `na`: time series of estimate roughness coefficient
- `x`: channel chainage
- `H`: time series of observed water surface elevation (reach)
- `W`: time series of observed width (reach)
- `S`: time series of observed slope (reach)
- `S0`: bed slope
- `hbf`: bankfull water surface elevation
- `wbf`: bankfull width
- `r`: channel shape parameter
- `z`: bed elevation

"""
function flow_parameters(Qa::Vector{Float64}, na::Vector{Float64}, x::Vector{Float64}, H::Vector{Float64}, W::Vector{Float64}, S::Vector{Float64}, S0::Vector{Float64}, hbf::Vector{Float64}, wbf::Vector{Float64}, r::Float64, z::Vector{Float64})
    ybf = hbf .- z
    ha = zeros(length(Qa))
    for t=1:length(Qa)
        h = gvf(Qa[t], (H[t]-z[1])*r/(r+1), S0, na[t], x, wbf, ybf, [r for _ in 1:length(x)])
        ha[t] = h[1]
    end
    A0 = (na .* Qa .* W.^(2/3) .* S.^(-1/2)).^(3/5) .- W .* ha
    mean(A0), mean(na)
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
- `nsamples`: sample size for rejection sampling

"""
function estimate(x::Vector{Float64}, H::Matrix{Float64}, W::Matrix{Float64}, Qp::Distribution, np::Distribution, rp::Distribution, zp::Distribution, nens::Int, nsamples::Int)
    wbf = maximum(W, dims=2)[:, 1]
    hbf = maximum(H, dims=2)[:, 1]
    S = diff(H, dims=1) ./ diff(x)
    S = [S[1, :]'; S]
    S0 = mean(S, dims=2)[:, 1]
    Qp, np, rp, zp = rejection_sampling(Qp, np, rp, zp, x, H, S0, wbf, hbf, nens, nsamples)
    Qe, ne, re, ze = prior_ensemble(x, Qp, np, rp, zp, nens)
    Se = repeat(S', outer=((nens ÷ size(H, 2)) + 1))'[:, 1:nens]
    bathymetry!(ze, Se, Qe, ne, re, x, hbf, wbf, mean(H, dims=2)[:, 1])
    zp = Truncated(Normal(mean(ze[1, :]), 1e-3), -Inf, minimum(H[1, :]))
    Sa = mean(Se, dims=2)[:, 1]
    Qa, Qu, na = assimilate(H, W, x, wbf, hbf, Sa, Qp, np, rp, zp, nens)
    za = zeros(length(x))
    za[1] = mean(zp)
    for i=2:length(x) za[i] = za[i+1] + Sa[i] * (x[i] - x[i-1]) end
    A0, n = flow_parameters(Qa, na, x, H[1, :], W[1, :], S[1, :], Sa, hbf, wbf, mean(re), za)
    Qa, Qu, A0, n
end

export
    Rectangular,
    Dingman,
    width,
    depth,
    gvf,
    rejection_sampling,
    prior_ensemble,
    assimilate


end
