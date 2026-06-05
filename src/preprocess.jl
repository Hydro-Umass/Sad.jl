using DataInterpolations
using Statistics
using BSplineKit
using LinearAlgebra

"""
    FloatM

Type alias for `Union{Missing, Float64}`. Used throughout SAD for SWOT
observation matrices that may contain missing values.
"""
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
    drop_unobserved(x, H, W, S)

Remove SWOT nodes that have no valid observations in either H or W across
all time steps. Resets chainage to start from zero.
"""
function drop_unobserved(x :: Vector{Float64}, H :: Matrix{FloatM}, W :: Matrix{FloatM}, S :: Matrix{FloatM})
    valid = [j for j in 1:size(H, 1)
             if !all(ismissing, H[j, :]) && !all(ismissing, W[j, :])]
    if isempty(valid)
        return x, H, W, S
    end
    x_out = x[valid] .- minimum(x[valid])
    return x_out, H[valid, :], W[valid, :], S[valid, :]
end

"""
    preprocess(xobs, Hobs, Wobs, Sobs; dx, min_slope)

Preprocess SWOT observations into a `SWOTReach` data structure.

Interpolates observations onto a regular computational chainage, estimates
bed slope, bankfull geometry, and builds PCHIP interpolants. Optionally
smooths WSE profiles with penalized B-splines.

# Arguments
- `xobs`:      chainage of SWOT observation nodes [m]
- `Hobs`:      WSE observations (nodes × timesteps), may contain `missing`
- `Wobs`:      width observations (nodes × timesteps), may contain `missing`
- `Sobs`:      slope observations (nodes × timesteps), may contain `missing`;
                pass `nothing` if slope data is unavailable

# Keyword arguments
- `dx`:         computational grid spacing [m] (default 200.0)
- `min_slope`:  minimum bed slope floor [-] (default 1e-5)
- `nknots`:     number of B-spline knots for WSE smoothing (default 15)
- `lambda`:     B-spline smoothing parameter; `nothing` = GCV selection (default)
"""
function preprocess(xobs :: Vector{Float64}, Hobs :: Matrix{FloatM}, Wobs :: Matrix{FloatM},
                    Sobs :: Union{Nothing, Matrix{FloatM}} = nothing;
                    dx :: Float64 = 200.0, min_slope :: Float64 = 1e-5,
                    nknots :: Int = 15, lambda :: Union{Nothing, Float64} = nothing,
                    outlier_thresh :: Float64 = 1.0)
    # If Sobs is not provided, create a dummy one for drop_unobserved
    if isnothing(Sobs)
        nobs_orig = size(Hobs, 1)
        Sobs_dummy = convert(Matrix{FloatM}, fill(1e-4, nobs_orig, size(Hobs, 2)))
    else
        Sobs_dummy = Sobs
    end
    xobs, Hobs, Wobs, Sobs_dummy = drop_unobserved(xobs, Hobs, Wobs, Sobs_dummy)
    nobs, nt = size(Hobs)
    x = build_chainage(xobs, dx)
    nx = length(x)

    # Return empty reach for degenerate cases
    if nobs < 2 || all(ismissing, Hobs) || all(ismissing, Wobs)
        itp_zero = PCHIPInterpolation(zeros(3), zeros(3); extrapolation=ExtrapolationType.Linear)
        return SWOTReach(
            SWOTObs(xobs, Hobs, Wobs),
            x,
            zeros(1, nt), zeros(1, nt),
            falses(nt),
            itp_zero, itp_zero, itp_zero, itp_zero,
            0.0, 1, 0, nt,
        )
    end

    # Bed slope and cumulative bed elevation from raw slope averaging
    S0obs = estimate_bed_slope(Sobs_dummy, min_slope)
    S0 = interpolate_to_chainage(xobs, S0obs, x)
    S0 = max.(S0, min_slope)
    S0_itp  = PCHIPInterpolation(S0,  x; extrapolation=ExtrapolationType.Linear)
    z = cumsum([0.0; S0[1:end-1] .* diff(x)])
    z_itp   = PCHIPInterpolation(z,   x; extrapolation=ExtrapolationType.Linear)

    # Bankfull WSE and width (median across timesteps, PCHIP to chainage)
    _row_median(row) = begin
        vals = collect(skipmissing(row))
        isempty(vals) ? 0.0 : median(vals)
    end
    wbf = interpolate_to_chainage(xobs, _row_median.(eachrow(Wobs)), x)
    wbf_itp = PCHIPInterpolation(wbf, x; extrapolation=ExtrapolationType.Linear)
    hbf = interpolate_to_chainage(xobs, _row_median.(eachrow(Hobs)), x)
    hbf_itp = PCHIPInterpolation(hbf, x; extrapolation=ExtrapolationType.Linear)

    # WSE and width at computational nodes per timestep
    # Outlier rejection: fit monotonic B-spline to WSE profile each timestep,
    # flag observations > outlier_thresh from the spline.
    if outlier_thresh > 0.0
        Hobs, Wobs = reject_wse_outliers(xobs, Hobs, Wobs, outlier_thresh;
                                         nknots=nknots, lambda=lambda,
                                         min_slope=min_slope)
    end
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
        if length(good_h) < 2 || length(good_w) < 2
            continue
        end
        _interp(y, x) = length(x) >= 3 ?
            PCHIPInterpolation(y, x; extrapolation=ExtrapolationType.Linear) :
            LinearInterpolation(y, x; extrapolation=ExtrapolationType.Linear)
        itp_w = _interp(convert(Vector{Float64}, w[good_w]), xobs[good_w])
        wout = itp_w.(xobs)
        itp_h = _interp(convert(Vector{Float64}, h[good_h]), xobs[good_h])
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
    reject_wse_outliers(xobs, Hobs, Wobs, threshold; kwargs...) -> (Hobs, Wobs)

Reject WSE (and corresponding width) observations that deviate more than
`threshold` meters from a monotonic smoothing B-spline fitted to the WSE
profile at each timestep.

For each timestep with valid WSE observations:
1. Fit a monotonic penalized B-spline to the WSE profile
2. Flag nodes where |WSE_obs - WSE_spline| > threshold
3. Set flagged WSE (and corresponding width) to `missing`

Thres prevents individual outlier nodes (vegetation, islands, instrument
errors) from biasing the Manning inference. The spline is only used for
outlier detection — spline values never enter inference.

# Arguments
- `xobs`: observation node locations [m]
- `Hobs`: WSE observations (nobs x nt) — modified in-place copy
- `Wobs`: width observations (nobs x nt) — modified in-place copy
- `threshold`: outlier threshold [m]; 0.0 disables rejection
- `nknots`: B-spline knots (default: max(8, nobs÷4))
- `lambda`: smoothing parameter; `nothing` = GCV (default)
- `min_slope`: minimum WSE slope for monotonicity enforcement

# Returns
(Hobs_clean, Wobs_clean) — copies with outliers set to `missing`
"""
function reject_wse_outliers(xobs::Vector{Float64},
                             Hobs::Matrix{FloatM},
                             Wobs::Matrix{FloatM},
                             threshold::Float64;
                             nknots::Int = max(8, length(xobs) ÷ 4),
                             lambda::Union{Nothing, Float64} = nothing,
                             min_slope::Float64 = 1e-5)
    nn, nt = size(Hobs)
    Hout = copy(Hobs)
    Wout = copy(Wobs)
    n_rejected = 0
    n_timesteps = 0

    for t in 1:nt
        # Get valid WSE obs at this timestep
        valid_h = findall(.!ismissing.(Hobs[:, t]))
        length(valid_h) < 5 && continue  # need enough nodes to fit a spline

        wse_vals = Float64.(Hobs[valid_h, t])
        x_vals = xobs[valid_h]

        # Fit monotonic smoothing spline
        # WSE must increase upstream (xobs sorted downstream→upstream)
        # In SAD convention: xobs[1] = downstream (lowest WSE), xobs[end] = upstream
        # So WSE should be monotonically increasing with xobs
        try
            spl = smooth_profile_bspline(x_vals, wse_vals;
                                         nknots=min(nknots, length(valid_h) ÷ 2),
                                         lambda=lambda,
                                         monotone_slope=min_slope)

            # Evaluate spline at observation nodes
            wse_spline = spl.(x_vals)

            # Flag outliers
            for (idx, node_idx) in enumerate(valid_h)
                deviation = abs(wse_vals[idx] - wse_spline[idx])
                if deviation > threshold
                    Hout[node_idx, t] = missing
                    Wout[node_idx, t] = missing
                    n_rejected += 1
                end
            end
            n_timesteps += 1
        catch
            # Spline fitting can fail for very sparse or noisy data
            # Just keep all observations in this case
            continue
        end
    end

    if n_rejected > 0
        total_obs = count(.!ismissing.(Hobs))
        @info "reject_wse_outliers: removed $n_rejected observations ($(round(n_rejected/total_obs*100, digits=1))%) across $n_timesteps timesteps (threshold=$(threshold)m)"
    end

    return Hout, Wout
end

"""
    interpolate_to_chainage(x, y, xnew)

Interpolate a complete vector to new chainage using PCHIP when 3+ points
are available, falling back to linear interpolation for 2 points.
"""
function interpolate_to_chainage(x :: Vector{Float64}, y :: Vector{Float64}, xnew :: Vector{Float64})
    if length(x) >= 3
        itp = PCHIPInterpolation(y, x; extrapolation=ExtrapolationType.Linear)
    else
        itp = LinearInterpolation(y, x; extrapolation=ExtrapolationType.Linear)
    end
    return itp.(xnew)
end

# ============================================================
# B-spline smoothing functions
# ============================================================

"""
    spline_design_matrix(x, B)

Build the B-spline design matrix A where A[i,j] = Bj(x_i).
"""
function spline_design_matrix(x::AbstractVector{<:Real}, B::BSplineBasis)
    n = length(x)
    m = length(B)
    A = zeros(n, m)
    for (i, xi) in enumerate(x)
        for (j, bj) in enumerate(B)
            A[i, j] = bj(xi)
        end
    end
    return A
end

"""
    penalty_matrix(n, order=2)

Build a penalty matrix for penalized B-spline fitting.
`order` specifies the order of the difference penalty (2 = second differences).
"""
function penalty_matrix(n::Int, order::Int=2)
    D = diff(I(n), dims=1)
    for _ in 2:order
        D = diff(D, dims=1)
    end
    return D
end

"""
    gcv_select(x, y, B; lambda_range)

Select the smoothing parameter λ that minimizes Generalized Cross-Validation.

GCV(λ) = RSS(λ) / (n - tr(H(λ)))² where RSS is the residual sum of squares
and H is the hat matrix. This balances fit quality against smoothness.

# Arguments
- `x`: observation locations
- `y`: observation values
- `B`: B-spline basis
- `lambda_range`: vector of λ values to search over

# Returns
The optimal λ value.
"""
function gcv_select(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                    B::BSplineBasis;
                    lambda_range::Vector{Float64} = [0.01, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0, 500.0])
    A = spline_design_matrix(x, B)
    D = penalty_matrix(length(B), 2)
    n = length(y)
    best_gcv = Inf
    best_lambda = lambda_range[1]
    for lam in lambda_range
        coeffs = (A'*A + lam*D'*D) \ (A'*y)
        yhat = A * coeffs
        resid = sum((y - yhat).^2)
        # Hat matrix trace
        H = A * ((A'*A + lam*D'*D) \ A')
        trH = tr(H)
        gcv = resid / (n - trH)^2
        if gcv < best_gcv
            best_gcv = gcv
            best_lambda = lam
        end
    end
    return best_lambda
end

"""
    smooth_profile_bspline(x, y; nknots, lambda, monotone_slope)

Fit a penalized cubic B-spline to scattered data with optional
monotonicity enforcement.

# Arguments
- `x`: observation locations (must be sorted in increasing order)
- `y`: observation values at each location
- `nknots`: number of internal knots for the B-spline basis (default: max(10, length(x)÷3))
- `lambda`: smoothing parameter (default: selected by GCV)
- `monotone_slope`: minimum slope to enforce (default: 0, no enforcement)

# Returns
A `Spline` object that can be evaluated at any x value.
"""
function smooth_profile_bspline(x::Vector{Float64}, y::Vector{Float64};
                                nknots::Int=max(10, length(x)÷3),
                                lambda::Union{Nothing, Float64}=nothing,
                                monotone_slope::Float64=0.0)
    n = length(x)
    n < 2 && error("smooth_profile_bspline: need at least 2 data points")

    # Create B-spline basis spanning the data range
    # Add boundary knots for proper interpolation at endpoints
    knots = collect(range(x[1], x[end], length=nknots))
    B = BSplineBasis(4, knots)  # Order 4 = cubic
    n_basis = length(B)
    n_basis > n && (B = BSplineBasis(4, collect(range(x[1], x[end], length=n÷2+2))); n_basis = length(B))

    # Design matrix
    A = spline_design_matrix(x, B)

    # Penalty matrix (2nd order differences)
    D = penalty_matrix(n_basis, 2)

    # Select smoothing parameter if not provided
    if isnothing(lambda)
        lambda = gcv_select(x, y, B)
    end

    # Fit penalized least squares
    coeffs = (A'*A + lambda*D'*D) \ (A'*y)
    spl = Spline(B, coeffs)

    # Enforce monotonicity if requested
    if monotone_slope > 0.0
        spl = enforce_monotonicity(spl, x, monotone_slope)
    end

    return spl
end

"""
    enforce_monotonicity(spl, x, min_slope)

Enforce that the B-spline derivative is at least `min_slope` at all
points in `x`. Adjusts the spline coefficients iteratively.

This is a simple post-processing step: evaluate the derivative at each
point, and if the slope is below `min_slope`, add a correction term
that ensures the minimum slope. Falls back to PCHIP if enforcement fails.
"""
function enforce_monotonicity(spl::Spline, x::Vector{Float64}, min_slope::Float64)
    B = basis(spl)
    dspl = BSplineKit.Splines._derivative(B, spl, Derivative(1))

    # Check if already monotone
    slopes = dspl.(x)
    if all(slopes .>= min_slope)
        return spl  # already monotone
    end

    # Enforce monotonicity by iteratively adjusting values
    # Simple approach: clip the values to be monotonically increasing
    # at minimum slope, then refit the spline
    y = spl.(x)
    for i in 2:length(x)
        required = y[i-1] + min_slope * (x[i] - x[i-1])
        if y[i] < required
            y[i] = required
        end
    end

    # Refit the spline to the adjusted values
    A = spline_design_matrix(x, B)
    D = penalty_matrix(length(B), 2)
    # Use a larger penalty to ensure smoothness of the corrected spline
    lambda = length(x) * 0.1  # moderate smoothing
    coeffs = (A'*A + lambda*D'*D) \ (A'*y)
    return Spline(B, coeffs)
end

"""
    slope_from_spline(spl, x)

Compute slope (first derivative) from a B-spline at the given locations.
Returns a vector of slope values with the same length as `x`.
"""
function slope_from_spline(spl::Spline, x::AbstractVector{<:Real})
    B = basis(spl)
    dspl = BSplineKit.Splines._derivative(B, spl, Derivative(1))
    slopes = dspl.(x)
    return max.(slopes, 1e-5)  # floor at minimum slope
end



