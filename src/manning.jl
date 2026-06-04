"""
    manning.jl — Analytical Manning + first-order GVF perturbation forward model

This module implements the SAD v2 forward model for water-surface elevation
prediction. It provides fully analytical, AD-compatible functions that replace
the numerical GVF ODE solver with a closed-form Manning solution augmented by
a first-order backwater perturbation.

## Forward Model: Manning + First-Order GVF Perturbation

Under the low-Froude-number regime (Fr² ≈ 0) typical of large low-gradient
rivers observed by SWOT, the GVF equation linearizes to an exponential backwater
decay from the downstream boundary condition. This yields a fully analytical
WSE profile with no ODE solve, making it both fast and AD-compatible.

### Formulation

1. **Uniform flow depth** (Manning + Dingman):
   y₀ = manning_depth(Q, n, r, Wb, Yb, S₀)

2. **Downstream boundary depth**:
   y_bc = (H_bc − z₀) · r/(r+1)

3. **Backwater length** (Dingman-corrected):
   λ = 3·y₀·r / ((10·r + 6)·S₀)
   Note: For r → ∞ (rectangular channel), λ → 3·y₀/(10·S₀)

4. **Depth profile**:
   y(x) = y₀ + (y_bc − y₀) · exp(−x/λ)

5. **WSE profile**:
   WSE(x) = z₀ + z_ref(x) + y(x) · (r+1)/r

### Key properties
- Fully analytical, no ODE solve, AD-compatible
- Captures dominant GVF effect (exponential backwater from downstream boundary)
- Smoothly interpolates: λ ≪ L → uniform flow; λ ≫ L → boundary-controlled
- Reduces to uniform-flow WSE when y_bc = y₀
"""

"""
    manning_Q(n, r, y, Wb, Yb, S) -> Float64

Analytical discharge from Manning's equation with the Dingman power-law
cross-section, under the wide-channel approximation (R ≈ y).

# Mathematical formulation
```julia
C = Wb * ((r+1) / (r*Yb))^(1/r)       # Dingman width-shape coefficient
Q = (1/n) * C * y^(5/3 + 1/r) * √S    # Manning discharge
```

Combining Dingman cross-sectional area
    A = Wb * ((r+1)/(r*Yb))^(1/r) * y^(1+1/r)
with Manning's equation Q = (1/n) * A * y^(2/3) * S^(1/2) yields the
combined exponent 5/3 + 1/r on depth.

# Arguments
- `n`:  Manning roughness coefficient [s/m^(1/3)]
- `r`:  Dingman channel shape exponent (1 < r < ∞)
- `y`:  mean flow depth [m]
- `Wb`: bankfull width [m]
- `Yb`: bankfull thalweg depth [m]
- `S`:  energy slope [-] (≈ bed slope for uniform flow)

# Returns
Discharge Q [m³/s]

# Examples
```jldoctest
julia> round(Sad.manning_Q(0.03, 2.0, 5.0, 200.0, 10.0, 1e-4), digits=1)
844.1
```
"""
function manning_Q(n::Real, r::Real, y::Real, Wb::Real, Yb::Real, S::Real)
    y = max(y, 0.01)  # depth positivity guard (see AGENTS.md conventions)
    C = Wb * ((r + 1) / (r * Yb))^(1 / r)
    return (1 / n) * C * y^(5 / 3 + 1 / r) * sqrt(S)
end

"""
    manning_depth(Q, n, r, Wb, Yb, S) -> Float64

Invert Manning's equation for mean flow depth under the Dingman cross-section
and wide-channel approximation.

# Mathematical formulation
```julia
C = Wb * ((r+1) / (r*Yb))^(1/r)       # Dingman width-shape coefficient
e = 5/3 + 1/r                         # depth exponent
y = (Q*n / (C * √S))^(1/e)            # Manning inverted for depth
```

This is the unique positive root of `manning_Q(n, r, y, Wb, Yb, S) - Q = 0`.

# Arguments
- `Q`:  discharge [m³/s]
- `n`:  Manning roughness coefficient [s/m^(1/3)]
- `r`:  Dingman channel shape exponent
- `Wb`: bankfull width [m]
- `Yb`: bankfull thalweg depth [m]
- `S`:  energy slope [-]

# Returns
Mean flow depth y [m]

# See also
[`manning_Q`](@ref) — the forward (Q from y) direction.
"""
function manning_depth(Q::Real, n::Real, r::Real, Wb::Real, Yb::Real, S::Real)
    Q = max(Q, 0.01)  # positivity guard
    C = Wb * ((r + 1) / (r * Yb))^(1 / r)
    e = 5 / 3 + 1 / r
    return (Q * n / (C * sqrt(S)))^(1 / e)
end

"""
    backwater_length(y0, S0, r) -> Float64

Compute the e-folding backwater length for the first-order GVF perturbation
with the Dingman cross-section.

# Mathematical formulation

Linearizing the GVF friction slope around uniform flow with Dingman geometry:
```julia
dSf/dy|_{y₀} = -(10/3 + 2/r) * S₀ / y₀
```

The exponential backwater decay length is:
```julia
λ = y₀ / ((10/3 + 2/r) * S₀)
  = 3*y₀*r / ((10*r + 6) * S₀)
```

For r → ∞ (rectangular channel with constant width), the 2/r term vanishes
and λ → 3·y₀/(10·S₀), the standard rectangular-channel backwater length.

# Physical interpretation
- λ ≪ L_reach: backwater effect is confined near the downstream boundary;
  most of the reach is in uniform flow.
- λ ≫ L_reach: backwater penetrates the entire reach; depth is controlled
  by the downstream boundary condition throughout.

# Arguments
- `y0`: uniform flow depth [m]
- `S0`: reach-averaged bed slope [-]
- `r`:  Dingman shape exponent

# Returns
Backwater length λ [m]
"""
function backwater_length(y0::Real, S0::Real, r::Real)
    y0 = max(y0, 0.01)
    return 3 * y0 * r / ((10 * r + 6) * S0)
end

"""
    backwater_depth(y0, y_bc, x, λ) -> Float64

Compute depth at chainage x under the first-order GVF perturbation.

    y(x) = y₀ + (y_bc − y₀) · exp(−x/λ)

At x = 0 (downstream boundary), y = y_bc. As x → ∞, y → y₀ (uniform flow).

# Arguments
- `y0`:  uniform flow depth [m]
- `y_bc`: downstream boundary depth [m]
- `x`:   chainage from downstream boundary [m]
- `λ`:   backwater length [m]

# Returns
Mean flow depth at chainage x [m]
"""
function backwater_depth(y0::Real, y_bc::Real, x::Real, λ::Real)
    return y0 + (y_bc - y0) * exp(-x / λ)
end

"""
    manning_wse(Q, n, r, z0, S, x_nodes, Wb_nodes, Yb_nodes, z_nodes)
        -> Vector

Predict the WSE profile at each computational node under uniform flow
(no backwater). Depth equals the reach-averaged Manning uniform-flow depth
everywhere.

Reach-averaged bankfull width and depth are computed as the mean of
`Wb_nodes` and `Yb_nodes` respectively. The slope `S` is assumed
reach-averaged.

# Mathematical formulation
```julia
Wb̄ = mean(Wb_nodes),  Ȳb = mean(Yb_nodes)
y₀  = manning_depth(Q, n, r, Wb̄, Ȳb, S)       # uniform flow depth
WSEⱼ = z₀ + z_refⱼ + y₀ * (r+1)/r            # WSE at each node
```

# Arguments
- `Q`:        discharge [m³/s]
- `n`:        Manning roughness coefficient
- `r`:        Dingman shape exponent
- `z0`:       downstream bed elevation [m]
- `S`:        reach-averaged bed slope [-]
- `x_nodes`:  computational chainage [m] (downstream=0, upstream=positive)
- `Wb_nodes`: bankfull width at each node [m]
- `Yb_nodes`: bankfull depth at each node [m]
- `z_nodes`:  cumulative bed elevation reference z_ref at each node [m]
              (z_ref[1] = 0 by convention; z₀ + z_ref(x) = total bed elevation)

# Returns
Vector of predicted WSE values at each node [m].
"""
function manning_wse(Q::Real, n::Real, r::Real, z0::Real, S::Real,
                     x_nodes::AbstractVector,
                     Wb_nodes::AbstractVector,
                     Yb_nodes::AbstractVector,
                     z_nodes::AbstractVector)
    Wb̄ = mean(Wb_nodes)
    Ȳb = mean(Yb_nodes)
    y0 = manning_depth(Q, n, r, Wb̄, Ȳb, S)
    y0 = max(y0, 0.01)
    return z0 .+ z_nodes .+ y0 * (r + 1) / r
end

"""
    manning_wse_backwater(Q, n, r, z0, S0, H_bc, x_nodes, Wb_nodes, Yb_nodes, z_nodes)
        -> Vector

Predict WSE profile including first-order GVF backwater perturbation.

Under the low-Froude-number assumption (Fr² ≈ 0), depth decays exponentially
from the downstream boundary depth toward the uniform-flow depth, with
e-folding length λ. This correctly captures M1 (backwater) and M2 (drawdown)
profiles.

# Mathematical formulation
```julia
Wb̄ = mean(Wb_nodes),  Ȳb = mean(Yb_nodes)
y₀  = manning_depth(Q, n, r, Wb̄, Ȳb, S₀)       # uniform flow depth
y_bc = (H_bc − z₀) * r/(r+1)                    # downstream boundary depth
λ   = backwater_length(y₀, S₀, r)               # e-folding length
y(x) = y₀ + (y_bc − y₀) · exp(−x/λ)             # depth profile
WSE(xⱼ) = z₀ + z_ref(xⱼ) + y(xⱼ) * (r+1)/r   # WSE at each node
```

**Important convention**: `H_bc` is the downstream boundary WSE — an INPUT to
the forward model, NOT an observation. In the SAD v2 likelihood, node j=1
(the downstream boundary) is excluded; the likelihood applies at nodes j ≥ 2
only.

# Boundary conditions
- At x = 0: WSE = H_bc (downstream boundary condition is satisfied exactly)
- At x → ∞: WSE → z₀ + z_ref + y₀·(r+1)/r (uniform flow)

# Arguments
- `Q`:        discharge [m³/s]
- `n`:        Manning roughness coefficient
- `r`:        Dingman shape exponent
- `z0`:       downstream bed elevation [m]
- `S0`:       reach-averaged bed slope [-]
- `H_bc`:     downstream boundary WSE [m] (typically `reach.H[1, t]`)
- `x_nodes`:  computational chainage [m]
- `Wb_nodes`: bankfull width at each node [m]
- `Yb_nodes`: bankfull depth at each node [m]
- `z_nodes`:  cumulative bed elevation reference at each node [m]

# Returns
Vector of predicted WSE values at each node [m].

# See also
[`manning_wse`](@ref) — uniform-flow profile (no backwater).
[`backwater_length`](@ref) — e-folding backwater distance.
"""
function manning_wse_backwater(Q::Real, n::Real, r::Real, z0::Real,
                               S0::Real, H_bc::Real,
                               x_nodes::AbstractVector,
                               Wb_nodes::AbstractVector,
                               Yb_nodes::AbstractVector,
                               z_nodes::AbstractVector)
    Wb̄ = mean(Wb_nodes)
    Ȳb = mean(Yb_nodes)
    y0 = manning_depth(Q, n, r, Wb̄, Ȳb, S0)
    y0 = max(y0, 0.01)
    y_bc = max((H_bc - z0) * r / (r + 1), 0.01)
    λ = backwater_length(y0, S0, r)
    y = y0 .+ (y_bc - y0) .* exp.(-x_nodes ./ λ)
    return z0 .+ z_nodes .+ y .* ((r + 1) / r)
end

# ============================================================
# Dingman cross-section geometry utilities
# ============================================================

"""
    area(y, Wb, Yb, r) -> Float64

Flow area for a Dingman power-law cross-section.

    A = Wb * (Ym / Yb)^(1/r) * y

where `Ym = (r+1)/r * y` is the thalweg depth and `y` is the mean depth.
"""
function area(y::Real, Wb::Real, Yb::Real, r::Real)
    Ym = (r + 1) / r * y
    Wb * (Ym / Yb)^(1 / r) * y
end

"""
    compute_wse(y, z, r) -> Float64

Compute water-surface elevation from mean flow depth.

    WSE = y * (r+1)/r + z

where `z` is the bed elevation at the same location.
"""
function compute_wse(y::Real, z::Real, r::Real)
    y * (r + 1) / r + z
end

"""
    compute_width(y, Wb, Yb, r) -> Float64

Compute water-surface width from mean flow depth.

    W = Wb * (y / Yb)^(1/r)

where `Yb` is the bankfull mean depth.
"""
function compute_width(y::Real, Wb::Real, Yb::Real, r::Real)
    Wb * (y / Yb)^(1 / r)
end