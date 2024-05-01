module Sad

abstract type CrossSection end

# river classification
@enum River braided=1 sinuous=2 more_sinuous=3 straight=4

# type for handling missing data
const FloatM = Union{Missing, Float64}

include("crosssections.jl")
include("gvf.jl")
include("priors.jl")
include("kalman.jl")

"""
    assimilate(H, W, S, x, wbf, hbf, S0, Qp, np, rp, zp, nens, constrain)

Assimilate SWOT observations for river reach.

# Arguments

- `H`: time series of water surface elevation profiles
- `W`: time series of water surface width profiles
- `S`: time series of water surface slope profiles
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
function assimilate(H, W, S, x, wbf, hbf, S0::Vector{Float64}, Qp::Distribution, np::Distribution, rp::Distribution, zp::Distribution, nens::Int, Hr::Union{Vector{FloatM}, Nothing}, constrain::Bool=true)
    seed!(1)
    min_ensemble_size = 5
    nt = size(H, 2)
    Qa = zeros(FloatM, nt)
    Qu = zeros(FloatM, nt)
    na = zeros(FloatM, nt)
    Qe, ne, re, ze = prior_ensemble(x, Qp, np, rp, zp, nens)
    for t in 1:nt
        if any(.!ismissing.(H[:, t])) # check if there are any valid observations
            hbc = ismissing(H[1, t]) ? interpolate_hbc(x, H[:, t], S0) : H[1, t]
            he = gvf_ensemble!(hbc, S0, x, hbf, wbf, Qe, ne, re, ze)
            # find valid flow depths simulated
            i = findall(he[1, :] .> 0)
            if length(i) < min_ensemble_size
                Sf = [!ismissing(S[j, t]) & (S[j, t] > 0) ? S[j, t] : minimum(S[:, t][.!ismissing.(S[:, t]) .& (S[:, t] .> 0)]) for j=1:length(x)]
                Wmean = mean.(skipmissing.(eachrow(W)))
                he = [(Qe[e] .* ne[e]) ./ ((ismissing(W[j, t]) ? Wmean[j] : W[j, t]) .* Sf[j].^0.5).^(3/5) for j=1:length(x), e=1:nens]
            end
            i = findall(he[1, :] .> 0)
            h = ze .+ he .* ((re .+ 1) ./ re)'
            X = zeros(2, length(i))
            X[1, :] = Qe[i]
            X[2, :] = ne[i]
            j = findall(.!ismissing.(H[:, t]))
            if isnothing(Hr)
                XA = h[j, i]
                d = convert(Vector{Float64}, H[j, t])
            else
                XA = zeros(1, length(i))
                XA[:] = mean(h[j, i], dims=1)
                d = [Hr[t]]
            end
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
                ahg_constrain!(Aq, hbc, hbf, wbf, S0, x, ne[i], re[i], ze[:, i])
            end
            Qa[t] = mean(Aq) > 0 ? mean(Aq) : mean(Aq[Aq .> 0])
            Qu[t] = std(Aq)
            An = A[2, :]
            na[t] = mean(An) > 0 ? mean(An) : mean(An[An .> 0])
        else
            Qa[t] = missing
            Qu[t] = missing
            na[t] = missing
        end
    end
    Qa, Qu, na
end

"""
    bathymetry!(ze, S, S0, Qe, ne, re, x, hbf, wbf, H, ϵₒ)

Estimate channel bed slope and thalweg (elevation) by assimilating SWOT water surface elevation observations.

# Arguments

- `ze`: ensemble of bed elevations
- `S`: time series of water surface slope profiles
- `S0`: bed slope
- `Qe`: ensemble of discharge
- `ne`: ensemble of roughness coefficient
- `re`: ensemble of channel shape parameters
- `x`: channel chainage
- `hbf`: bankfull water surface elevation
- `wbf`: bankfull width
- `H`: time-averaged water surface elevation observed
- `ϵₒ`: observation error (default value of 10 cm)

"""
function bathymetry!(ze::Matrix{Float64}, S::Matrix{FloatM}, S0::Vector{Float64}, Qe::Vector{Float64}, ne::Vector{Float64}, re::Vector{Float64}, x::Vector{Float64}, hbf::Vector{Float64}, wbf::Vector{Float64}, H::Vector{Float64}, ϵₒ::Float64=0.1)
    seed!(1)
    nens = size(ze, 2)
    Se = repeat(S', outer=((nens ÷ size(H, 2)) + 1))'[:, 1:nens]
    for k=1:length(x)
        Se[k, ismissing.(Se[k, :])] .= S0[k]
        Se[k, Se[k, :] .< 0] .= S0[k]
    end
    Se = convert(Matrix{Float64}, Se)
    min_ensemble_size = 5
    he = gvf_ensemble!(H[1], Se, x, hbf, wbf, Qe, ne, re, ze)
    i = findall(he[1, :] .> 0)
    if length(i) > min_ensemble_size
        h = ze .+ he .* ((re .+ 1) ./ re)'
        X = ze[:, i]
        XA = h[:, i]
        d = H
        E = (((d .- mean(XA,dims=2)).^2 .- std(XA, dims=2).^2))
        E[E .> ϵₒ] .= ϵₒ
        A = letkf(X, d, XA, E, diagR=true)
        ze[:, i] .= A
        S = diff(ze, dims=1) ./ diff(x)
        S = [S[1, :]'; S]
        for k=1:length(x) S[k, S[k, :] .< 0] .= minimum(Se[k, :]) end
        Se = S
    end
    Se
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
    flow_parameters(Qa, na, W, S, dA)

Estimate flow parameters (roughness coefficient and baseflow cross-sectional area) from assimilated discharge.

# Arguments

- `Qa`: time series of estimated discharge
- `na`: time series of estimate roughness coefficient
- `W`: time series of observed width (reach)
- `S`: time series of observed slope (reach)
- `dA`: time series of changes in reach-averaged cross-sectional area (from minimum observed)

"""
function flow_parameters(Qa::Vector{FloatM}, na::Vector{FloatM}, Wm::Vector{FloatM}, Sm::Vector{FloatM}, dA::Vector{FloatM})
    i = findall(Sm .> 0)
    A0 = (na[i] .* Qa[i]).^(3/5) .* Wm[i].^(2/5) .* Sm[i].^(-3/10) .- dA[i]
    A0 = mean(skipmissing(A0))
    n = mean(skipmissing(na))
    # catch implausible negative value of A0
    A0 > 0 ? A0 : 0.0, n
end

"""
    drop_unobserved(x, H, W, S)

Remove cross sections with no valid observations.

# Arguments

- `x`: channel chainage
- `H`: time series of water surface elevation profiles
- `W`: time series of width profiles
- `S`: time series of water surface slope profiles

"""
function drop_unobserved(x::Vector{Float64}, H::Matrix{FloatM}, W::Matrix{FloatM}, S::Matrix{FloatM})
    i = [j for j=1:size(H, 1) if (!all(ismissing.(H[j, :])) & !all(ismissing.(W[j, :])))]
    if isempty(i)
        x, H, W, S
    else
        x[i] .- minimum(x[i]), H[i, :], W[i, :], S[i, :]
    end
end

"""
    calc_bed_slope(x, S)

Estimate initial channel bed slope.

# Arguments

- `x`: channel chainage
- `S`: time series of water surface slope profiles

"""
function calc_bed_slope(x::Vector{Float64}, S::Matrix{FloatM})
    # calculate initial bed slope as the mean water surface slope
    # when data are missing, we may end up with negative slopes which
    # can cause numerical instabilities in the GVF model so we use either
    # the minimum positive slope or the entire reach time-averaged slope
    S0 = zeros(length(x))
    Sm = mean.(skipmissing.(eachrow(S)))
    for k=1:length(x)
        s = S[k, :]
        s = s[.!ismissing.(s)]
        if length(s) < 1 || maximum(s) < 0
            S0[k] = mean(Sm[.!isnan.(Sm)])
        else
            S0[k] = mean(s[s .> 0])
        end
    end
    S0
end

"""
    calc_dA(Qa, na, W, S)

Calculate changes in cross-sectional area.

# Arguments

- `Qa`: time series of estimated discharge
- `na`: time series of estimated roughness coefficient
- `W`: time series of width profiles
- `S`: time series of water surface slope profiles

"""
function calc_dA(Qa::Vector{FloatM}, na::Vector{FloatM}, W::Matrix{FloatM}, S::Matrix{FloatM})
    Wm = mean.(skipmissing.(eachcol(W)))
    Sm = mean.(skipmissing.(eachcol(S)))
    Sm[Sm .< 0] .= NaN
    A = (na .* Qa .* Wm.^(2/3) .* Sm.^(-1/2)).^(3/5)
    A = [!ismissing(a) & isnan(a) ? missing : a for a in A]
    Amin = minimum(skipmissing(A))
    dA = convert(Vector{FloatM}, A .- Amin)
    dA
end

"""
    estimate(x, H, W, S, Qp, np, rp, zp, nens)

Estimate discharge and flow parameters from SWOT observations.

# Arguments
- `x`: channel chainage
- `H`: time series of water surface elevation profiles
- `W`:time series of width profiles
- `S`: time series of water surface slope profiles
- `Qp`: discharge prior distribution
- `np`: roughness coefficient prior distribution
- `rp`: channel shape parameter prior distribution
- `zp`: downstream bed elevation prior distribution
- `nens`: ensemble size
- `nsamples`: sample size for rejection sampling

"""
function estimate(x::Vector{Float64}, H::Matrix{FloatM}, W::Matrix{FloatM}, S::Matrix{FloatM}, dA::Union{Vector{FloatM}, Nothing}, Qp::Distribution, np::Distribution, rp::Distribution, zp::Distribution, nens::Int, nsamples::Int, Hr::Union{Vector{FloatM}, Nothing}=nothing, Wr::Union{Vector{FloatM}, Nothing}=nothing, Sr::Union{Vector{FloatM}, Nothing}=nothing)
    wbf = maximum.(skipmissing.(eachrow(W)))
    hbf = maximum.(skipmissing.(eachrow(H)))
    S0 = calc_bed_slope(x, S)
    Qp, np, rp, zp = rejection_sampling(Qp, np, rp, zp, x, H, S0, wbf, hbf, nens, nsamples)
    Qe, ne, re, ze = prior_ensemble(x, Qp, np, rp, zp, nens)
    Se = bathymetry!(ze, S, S0, Qe, ne, re, x, hbf, wbf, mean.(skipmissing.(eachrow(H))))
    zp = truncated(Normal(mean(ze[1, :]), 1e-3), -Inf, minimum(skipmissing(H[1, :])))
    Sa = mean(Se, dims=2)[:, 1]
    Qa, Qu, na = assimilate(H, W, S, x, wbf, hbf, Sa, Qp, np, rp, zp, nens, Hr, false)
    # diagnose discharge and roughness coefficient
    Qa = convert(Vector{FloatM}, [clamp(q, minimum(Qp), maximum(Qp)) for q in Qa])
    nam = [n for n in na if !ismissing(n) & !isnan(n)]
    nm = length(nam) > 0 ? mean(nam) : mean(ne)
    na = convert(Vector{FloatM}, [!ismissing(n) & isnan(n) ? nm : n for n in na])
    za = zeros(length(x))
    za[1] = mean(zp)
    for i=2:length(x) za[i] = za[i-1] + Sa[i] * (x[i] - x[i-1]) end
    if isnothing(dA)
        dA = calc_dA(Qa, na, W, S)
    end
    Wm = isnothing(Wr) ? mean.(skipmissing.(eachcol(W))) : Wr
    Sm = isnothing(Sr) ? mean.(skipmissing.(eachcol(S))) : Sr
    A0, n = flow_parameters(Qa, na, convert(Vector{FloatM}, Wm), convert(Vector{FloatM}, Sm), dA)
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
