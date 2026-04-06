using DataInterpolations
using Statistics

# type alias for float or missing data
const FloatM = Union{Missing, Float64}

"""
SWOT observations data structure.
- `x`: chainage of SWOT nodes (nodes)
- `H`: WSE observations (nodes x times)
- `W`: width observations (nodes x times)
"""
struct SWOTObs
    x :: Vector{Float64}
    H :: Matrix{FloatM}
    W :: Matrix{FloatM}
end

"""
SWOT reach data structure
- `obs`: SWOT observations
- `x`: computational chainage, downstream → upstream
- `wbf`: interpolant for bankfull width
- `hbf`: interpolant for bankfull WSE
- `S0`: interpolant for bed slope
- `z`: interpolant for bed elevation
- `H`: monotonic WSE interpolated to `x`
- `W`: width interpolated to `x`
- `valid`: bit vector with true where timestep has ≥ 2 valid observations
- `hmin`: minimum observed WSE at downstream node
- `nx`: number of computational nodes
- `nobs`: number of SWOT nodes with at least one valid observation
- `nt`: number of timesteps
"""
struct SWOTReach{I}
    obs :: SWOTObs
    x :: Vector{Float64}
    H :: Matrix{Float64}
    W :: Matrix{Float64}
    valid :: BitVector
    S0 :: I
    wbf :: I
    hbf :: I
    z :: I
    hmin :: Float64
    nx :: Int
    nobs :: Int
    nt :: Int
end

"""
    drop_unobserved!(x, H, W, S)

Remove SWOT nodes that have no valid observations in either H or W across
all time steps. Resets chainage to start from zero.
"""
function drop_unobserved!(x :: Vector{Float64}, H :: Matrix{FloatM}, W :: Matrix{FloatM}, S :: Matrix{FloatM})
    valid = [j for j in 1:size(H, 1)
             if !all(ismissing, H[j, :]) && !all(ismissing, W[j, :])]
    if isempty(valid)
        return x, H, W, S
    end
    x_out = x[valid] .- minimum(x[valid])
    return x_out, H[valid, :], W[valid, :], S[valid, :]
end

"""
        preprocess(xobs, Hobs, Wobs, Sobs, dx, min_slope)
"""
function preprocess(xobs :: Vector{Float64}, Hobs :: Matrix{FloatM}, Wobs :: Matrix{FloatM}, Sobs :: Matrix{FloatM}; dx :: Float64 = 200.0, min_slope :: Float64 = 1e-5)
    xobs, Hobs, Wobs, Sobs = drop_unobserved!(xobs, Hobs, Wobs, Sobs)
    nobs, nt = size(Hobs)
    # if all SWOT observations are missing after dropping unobserved nodes,
    # return an empty reach with no valid timesteps
    if nobs == 0 || all(ismissing, Hobs) || all(ismissing, Wobs)
        x_empty  = [0.0]
        itp_zero = PCHIPInterpolation([0.0], [0.0]; extrapolation=ExtrapolationType.Linear)
        return SWOTReach(
            SWOTObs(xobs, Hobs, Wobs),
            x_empty,
            zeros(1, nt), zeros(1, nt),
            falses(nt),
            itp_zero, itp_zero, itp_zero, itp_zero,
            0.0, 1, 0, nt,
        )
    end
    x = build_chainage(xobs, dx)
    nx = length(x)
    S0obs = estimate_bed_slope(Sobs, min_slope)
    S0 = interpolate_to_chainage(xobs, S0obs, x)
    S0 = max.(S0, min_slope)
    S0_itp  = PCHIPInterpolation(S0,  x; extrapolation=ExtrapolationType.Linear)
    z = cumsum([0.0; S0[1:end-1] .* diff(x)])
    z_itp   = PCHIPInterpolation(z,   x; extrapolation=ExtrapolationType.Linear)
    # use median across timesteps for bankfull width and WSE as maximum causes
    # GVF to produce near-zero depths
    # falls back to zero if no valid SWOT observations were found
    _row_median(row) = begin
        vals = collect(skipmissing(row))
        isempty(vals) ? 0.0 : median(vals)
    end
    wbf = interpolate_to_chainage(xobs, _row_median.(eachrow(Wobs)), x)
    wbf_itp = PCHIPInterpolation(wbf, x; extrapolation=ExtrapolationType.Linear)
    hbf = interpolate_to_chainage(xobs, _row_median.(eachrow(Hobs)), x)
    hbf_itp = PCHIPInterpolation(hbf, x; extrapolation=ExtrapolationType.Linear)
    H_f, W_f, valid = obs_chainage(xobs, Hobs, Wobs, min_slope)
    W = zeros(nx, nt)
    H = zeros(nx, nt)
    for t in findall(valid)
        H[:, t] = interpolate_to_chainage(xobs, H_f[:, t], x)
        W[:, t] = interpolate_to_chainage(xobs, W_f[:, t], x)
    end
    hmin = minimum(skipmissing(Hobs[1, :]))
    SWOTReach(SWOTObs(xobs, Hobs, Wobs), x, H, W, valid, S0_itp, wbf_itp, hbf_itp, z_itp, hmin, nx, nobs, nt)
end

"""
    build_chainage(xobs, dx)

Build a uniform chainage for the river reach. Starts downstream (x=0)
and ends at the most upstream SWOT node.
"""
function build_chainage(xobs :: Vector{Float64}, dx :: Float64)
    xmax = maximum(xobs)
    n = ceil(Int, xmax / dx) + 1
    return range(0.0, step=dx, length=n) |> collect
end

"""
    estimate_bed_slope(Sobs :: Matrix{Float64}, min_slope :: Float64)

Estimate the bed slope profile of the river reach from the time-averaged SWOT water surface slope.
"""
function estimate_bed_slope(Sobs :: Matrix{FloatM}, min_slope :: Float64)
    n = size(Sobs, 1)
    S0 = zeros(n)
    all_valid = [s for s in skipmissing(vec(Sobs)) if s > 0]
    S_reach = isempty(all_valid) ? min_slope : mean(all_valid)
    for k in 1:n
        row = Sobs[k, :]
        valid  = [s for s in skipmissing(row) if s > 0]
        S0[k] = isempty(valid) ? S_reach : mean(valid)
    end
    return S0
end

"""
    obs_chainage(xobs, Hobs, Wobs, x, min_slope)

Fill gaps in observed WSE and width profiles via PCHIP interpolation and enforce monotonicity in WSE.
 Return `nothing` if fewer than two valid observations are available.
"""
function obs_chainage(xobs :: Vector{Float64}, Hobs :: Matrix{FloatM}, Wobs :: Matrix{FloatM}, min_slope :: Float64)
    nn, nt = size(Hobs)
    Hout = zeros(nn, nt)
    Wout = zeros(nn, nt)
    valid = falses(nt)
    for t in 1:nt
        h = Hobs[: ,t]
        w = Wobs[:, t]
        good_h = findall(.!ismissing.(h))
        good_w = findall(.!ismissing.(w))
        if length(good_h) < 2
            continue
        end
        itp_w = PCHIPInterpolation(convert(Vector{Float64}, w[good_w]), xobs[good_w]; extrapolation=ExtrapolationType.Linear)
        wout = itp_w.(xobs)
        itp_h = PCHIPInterpolation(convert(Vector{Float64}, h[good_h]), xobs[good_h]; extrapolation=ExtrapolationType.Linear)
        hout = itp_h.(xobs)
        for k in 2:nn
            dx = xobs[k] - xobs[k-1]
            required = hout[k-1] + min_slope * dx
            if hout[k] < required
                hout[k] = required
            end
        end
        Hout[:, t] = hout
        Wout[:, t] = wout
        valid[t] = true
    end
    return Hout, Wout, valid
end

"""
    interpolate_to_chainage(x, y, xnew)

Interpolate a complete vector to new chainage using PCHIP.
"""
function interpolate_to_chainage(x :: Vector{Float64}, y :: Vector{Float64}, xnew :: Vector{Float64})
    itp = PCHIPInterpolation(y, x; extrapolation=ExtrapolationType.Linear)
    return itp.(xnew)
end
