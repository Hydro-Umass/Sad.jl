using DifferentialEquations
using DataInterpolations

# numerical regularization term for the GVF denominator (1 - Fr² + ε)
# prevents numerical instability if the solver gets near critical flow.
# at the reach scales of SWOT (large, low-gradient rivers) Fr << 1,
# so this constant has no physical effect on results
const gvf_ε = 1e-6

"""
Flow area for a Dingman cross section.

    A = Wb * (Ym / Yb)^(1/r) * y

where y is the mean flow depth and Ym = (r+1)/r * y is the thalweg depth.
"""
function area(y::Real, Wb::Real, Yb::Real, r::Real)
    Ym = (r + 1) / r * y
    Wb * (Ym / Yb)^(1 / r) * y
end

"""
    gvf_rhs(y, p, x)

ODE right-hand side for the Gradually-Varied-Flow equation, integrated
upstream (x increases from downstream to upstream boundary).

    dy/dx = -(S0(x)·α - Sf(y)) / (1 - Fr(y)² + ε)

The negative sign arises because we integrate in the upstream direction
(increasing x): WSE rises upstream so dy/dx > 0 for normal backwater, which
requires S0 > Sf. The standard downstream-integration form has the opposite
sign convention.

Manning friction slope:  Sf = (n·Q)² / (A²·y^(4/3))
Froude number:           Fr = Q / (A·√(g·y))

Both use the wide-channel approximation R ≈ y (mean depth), standard for
large low-gradient rivers observed by SWOT.

# Parameter tuple `p`
- `p[1]` Q:       discharge [m³/s]
- `p[2]` n:       Manning roughness coefficient
- `p[3]` r:       Dingman shape exponent
- `p[4]` α:       slope correction factor (S_actual = S0·α)
- `p[5]` S0_itp:  interpolant for bed slope S0(x)
- `p[6]` wbf_itp: interpolant for bankfull width Wb(x) [m]
- `p[7]` hbf_itp: interpolant for bankfull WSE hbf(x) [m]
- `p[8]` z_itp:   interpolant for cumulative bed elevation reference z(x) [m]
- `p[9]` z0:      downstream bed elevation [m]
"""
function gvf_rhs(y, p, x)
    Q, n, r, α, S0_itp, wbf_itp, hbf_itp, z_itp, z0 = p
    y  = max(y, 0.01)
    Wb = wbf_itp(x)
    # bed elevation at x: reference profile scaled by α, anchored at z0
    # z(x) = z0 + z_ref(x)·α,  where z_ref[1] = 0 by construction
    zx = z0 + z_itp(x) * α
    Yb = max(hbf_itp(x) - zx, 0.01)
    S0 = S0_itp(x) * α
    A  = area(y, Wb, Yb, r)
    Sf = (n * Q / (A * y^(2/3)))^2
    Fr = Q / (A * sqrt(9.806 * y))
    -(S0 - Sf) / (1 - Fr^2 + gvf_ε)
end

"""
    gvf_solve(Q, n, r, α, z0, H_bc, reach; saveat) → Vector{Float64} or nothing

Solve the GVF equation for a single parameter set and return predicted
water surface elevation at the locations specified by `saveat`.

# Arguments
- `Q`:      discharge [m³/s]
- `n`:      Manning roughness coefficient
- `r`:      Dingman channel shape exponent
- `α`:      slope correction factor (S_actual = S0·α)
- `z0`:     downstream bed elevation [m]
- `H_bc`:   downstream boundary WSE [m]  (use `reach.H[1, t]`)
- `reach`:  `SWOTReach` from preprocessing stage
- `saveat`: chainage locations at which to return WSE [m]
            (default: `reach.obs.x`, aligning output with SWOT observations)

# Returns
Predicted WSE vector at `saveat` locations, or `nothing` if the solve fails.
"""
function gvf_solve(Q::Real, n::Real, r::Real, α::Real, z0::Real, H_bc::Real,
                   reach::SWOTReach;
                   saveat::Vector{Float64} = collect(reach.obs.x))
    x       = reach.x
    S0_itp  = reach.S0
    wbf_itp = reach.wbf
    hbf_itp = reach.hbf
    z_itp   = reach.z
    y_bc    = max((H_bc - z0) * r / (r + 1), 0.01)
    p       = (Q, n, r, α, S0_itp, wbf_itp, hbf_itp, z_itp, z0)
    prob    = ODEProblem(gvf_rhs, y_bc, (x[1], x[end]), p)
    sol     = try
        solve(prob, Tsit5(),
              abstol=1e-3, reltol=1e-3,
              saveat=saveat, maxiters=10_000)
    catch
        return nothing
    end
    sol.retcode == ReturnCode.Success || return nothing
    z_saveat = z0 .+ z_itp.(saveat) .* α
    sol.u .* ((r + 1) / r) .+ z_saveat
end
