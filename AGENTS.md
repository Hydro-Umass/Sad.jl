# SAD.jl — Agent Context

## Project Overview

SAD (SWOT Assimilated Discharge) estimates river discharge from SWOT satellite observations (WSE, width, slope). The package implements a data assimilation algorithm that combines an analytical hydraulic forward model with prior distributions to infer discharge at each river reach and timestep.

**Current version (v2)** is the active implementation. The v1 two-stage LETKF+GVF pipeline has been fully replaced.

---

## SAD v2 Architecture

### Forward Model: Manning + First-Order GVF Perturbation

Under the low-Froude-number regime (Fr² ≈ 0), the GVF equation linearizes to an exponential backwater decay:

```
y₀  = manning_depth(Q, n, r, S₀, Wb, Yb)           # uniform flow depth, analytical
y_bc = (H_bc − z₀) · r/(r+1)                         # downstream boundary depth
λ   = 3·y₀·r / ((10·r + 6)·S₀)                     # backwater length (Dingman-corrected)
y(x) = y₀ + (y_bc − y₀) · exp(−x/λ)                 # depth at node x, analytical
WSE(x) = z₀ + z_ref(x) + y(x) · (r+1)/r             # WSE prediction at node x
```

Key properties:
- Fully analytical, no ODE solve, AD-compatible (ForwardDiff throughout)
- Captures dominant GVF effect (exponential backwater decay from downstream boundary)
- Smoothly interpolates: λ ≪ L_reach → uniform flow; λ ≫ L_reach → boundary-controlled
- Reduces to uniform-flow WSE when y_bc = y₀

### Inference: Joint MAP + Laplace Uncertainty

All parameters estimated simultaneously by minimizing a single objective:

- **Per-node WSE likelihood**: Student-t (ν ≈ 3–5) for robustness to outliers. Applied at upstream nodes only — the downstream WSE (H_bc) is a boundary condition, not a data point
- **Per-node width likelihood** (optional): Gaussian likelihood on predicted vs observed width, activated via `use_width=true`. Uses `σ_W` estimated from data
- **Q priors**: monthly climatological log-normal priors from SoS database (σ=1.0 in log-space, giving wide 95% CI ≈ [0.14, 7.4]×median)
- **Static parameter priors**: truncated Normal on n ∈ [0.01, 0.05], truncated Normal/Uniform on r, **truncated Normal on z₀** centred at Manning depth estimate with upper bound below hmin
- **Temporal smoothness**: penalty on |log Q_t − log Q_{t-1}|², adaptively scaled for data-sparse reaches
- **Optimizer**: L-BFGS with ForwardDiff autodiff (~40+3 parameters, converges in ~50–200 iterations)
- **Uncertainty**: Laplace approximation — inverse Hessian at MAP gives posterior covariance, Q_std via delta method

### Data-Adaptive σ_obs and Retry Logic

WSE observation noise σ_obs is auto-estimated from local residuals in the observations (with a 0.3 m floor). If the optimizer converges to unphysical parameters (n > 0.10, r > 15, z₀ far from hmin, or Q spikes > 10× prior median), σ_obs is progressively doubled until convergence. If after retry the final n exceeds the prior truncation (0.05), a prior+WSE anomaly fallback is used instead — this handles cases where the Manning+Dingman model has structural error too large for the estimated σ_obs.

### Prior + WSE Anomaly Fallback

When MAP inference fails to produce physical parameters even after σ_obs retry, the fallback computes:
- `Q_t = Q_prior(t) × α_t` where α_t is a damped WSE depth anomaly (β=1.5, blended 50% toward prior)
- Static parameters set to prior medians
- Uncertainty set to prior standard deviation
- Flagged with `fallback=true` in the result

---

## Key Data Structures

### `SWOTObs`
Defined in `src/preprocess.jl`. Raw SWOT observations:
- `x::Vector{Float64}` — chainage of SWOT observation nodes
- `H::Matrix{FloatM}` — WSE observations (nodes × timesteps with missings)
- `W::Matrix{FloatM}` — width observations (nodes × timesteps with missings)

### `SWOTReach`
Defined in `src/preprocess.jl`. Preprocessed reach data:
- `obs::SWOTObs` — raw SWOT observations
- `x::Vector{Float64}` — computational chainage (downstream to upstream, x[1]=0)
- `H::Matrix{Float64}`, `W::Matrix{Float64}` — interpolated WSE/width on computational grid
- `valid::BitVector` — timesteps with ≥2 valid observations
- `S0, wbf, hbf, z::I` — PCHIP interpolants for bed slope, bankfull width/WSE, bed elevation
- `hmin::Float64` — minimum observed downstream WSE
- `nx, nobs, nt::Int` — grid/node/timestep counts

### `SWOTPriors`
Defined in `src/priors.jl`. Prior distributions:
- `Qp` — discharge prior (single `Distribution` or 12-element `Vector{Distribution}` for monthly)
- `np::Distribution` — Manning n: truncated Normal(0.03, 0.01, [0.01, 0.05])
- `rp::Distribution` — Dingman r: truncated Normal or Uniform depending on data source
- `zp::Distribution` — bed elevation z₀: truncated Normal centred at Manning depth estimate, upper bound below hmin (minimum depth 0.5 m floor)

### `ManningPrecomp`
Defined in `src/inference.jl`. Pre-computed reach data for the Manning forward model:
- `x_nodes::Vector{Float64}` — SWOT observation node chainage
- `Wb_nodes::Vector{Float64}` — bankfull width at each observation node
- `hbf_z_nodes::Vector{Float64}` — hbf(x) − z_ref(x) at each node (z₀-independent part of Yb)
- `z_ref_nodes::Vector{Float64}` — cumulative bed elevation reference z_ref
- `S0::Float64` — reach-averaged bed slope
- `nt::Int` — total number of timesteps
- `valid_ts::Vector{Int}` — indices of valid timesteps
- `H_bc::Vector{Float64}` — downstream boundary WSE per valid timestep
- `upstream_j::Vector{Vector{Int}}` — non-missing upstream node indices per timestep
- `upstream_H::Vector{Vector{Float64}}` — observed WSE at upstream nodes
- `upstream_W::Vector{Vector{Float64}}` — observed width at upstream nodes
- `months_v::Vector{Int}` — calendar month per valid timestep
- `σ_obs_est::Float64` — auto-estimated WSE noise
- `σ_W_est::Float64` — auto-estimated width noise (fractional uncertainty)

---

## Key Physics

### Dingman Cross-Section
Power-law channel shape parameterized by r (1 < r < ∞):
- `A = Wb · (Ym/Yb)^(1/r) · y` — cross-sectional area
- `W = Wb · (y/Yb)^(1/r)` — water surface width
- `Y = y · (r+1)/r + z` — WSE from mean depth (and inverse: `y = (WSE - z) · r/(r+1)`)
- `Ym = (r+1)/r · y` — thalweg depth; `Yb = (r+1)/r · yb` — bankfull thalweg depth

### Manning's Equation (wide-channel, R ≈ y)
```
Q = (1/n) · A · y^(2/3) · S^(1/2)
```
With Dingman area: `A = Wb · ((r+1)/(r·Yb))^(1/r) · y^(1+1/r)`, so:
```
Q = (1/n) · Wb · ((r+1)/(r·Yb))^(1/r) · y^(5/3+1/r) · S^(1/2)
```
Inverted for depth: `y = (Q·n / (Wb·((r+1)/(r·Yb))^(1/r)·√S))^(1/(5/3+1/r))`

### Backwater Length (Dingman-corrected)
```
λ = 3·y₀·r / ((10·r + 6)·S₀)
```
For r → ∞ (rectangular channel), λ → 3·y₀/(10·S₀). The Dingman correction factor 2/r in the friction slope derivative reduces the backwater length for finite r.

### First-Order GVF Perturbation
```
y(x) = y₀ + (y_bc − y₀) · exp(−x/λ)
```
At x=0 (downstream boundary), y = y_bc. As x → ∞, y → y₀ (uniform flow).

---

## File Structure

```
src/
  Sad.jl              # Module — includes all source files
  preprocess.jl       # SWOT data preprocessing, SWOTObs/SWOTReach structs,
                      #   PCHIP interpolation (DataInterpolations.jl), B-spline
                      #   smoothing, bed slope estimation, chainage construction
  priors.jl           # SWOTPriors struct, SoS prior construction, monthly Q priors,
                      #   z₀ prior from Manning depth estimate, uninformative priors
  manning.jl          # Analytical Manning + first-order GVF perturbation forward model:
                      #   manning_Q, manning_depth, backwater_length, backwater_depth,
                      #   manning_wse, manning_wse_backwater, plus Dingman geometry
                      #   utilities (area, compute_wse, compute_width)
  inference.jl        # Joint MAP inference:
                      #   ManningPrecomp, precompute_manning, estimate_σ_obs,
                      #   estimate_σ_W, neg_log_joint, infer, _prior_anomaly_fallback,
                      #   initialize_theta, laplace_uncertainty, safe_logpdf
```

---

## Key Functions

### `manning.jl`
| Function | Description |
|---|---|
| `manning_Q(n, r, y, Wb, Yb, S)` | Analytical discharge from Manning+Dingman |
| `manning_depth(Q, n, r, Wb, Yb, S)` | Invert Manning for mean depth |
| `backwater_length(y0, S0, r)` | E-folding backwater decay distance |
| `backwater_depth(y0, y_bc, x, λ)` | Depth at chainage x under GVF perturbation |
| `manning_wse(Q, n, r, z0, S, x_nodes, Wb_nodes, Yb_nodes, z_nodes)` | WSE profile (uniform flow) |
| `manning_wse_backwater(Q, n, r, z0, S0, H_bc, x_nodes, Wb_nodes, Yb_nodes, z_nodes)` | WSE profile with backwater |
| `area(y, Wb, Yb, r)` | Dingman cross-sectional area |
| `compute_wse(y, z, r)` | WSE from mean depth |
| `compute_width(y, Wb, Yb, r)` | Water surface width from mean depth |

### `inference.jl`
| Function | Description |
|---|---|
| `precompute_manning(reach, months; S0)` | Extract geometry/observations for forward model |
| `estimate_σ_obs(x_nodes, upstream_j, upstream_H, S0)` | Auto-estimate WSE noise from observations |
| `estimate_σ_W(upstream_W)` | Auto-estimate width noise from observations |
| `neg_log_joint(θ, precomp, priors; σ_obs, ν, λ_smooth, use_width)` | MAP objective function |
| `infer(priors, reach; time_str, σ_obs, ν, λ_smooth, ...)` | Main inference entry point: joint MAP, σ_obs retry, profile refit, Laplace uncertainty |
| `_profile_refit(θ_map, precomp, priors; ...)` | Fix n, r at MAP and re-optimize Q, z₀ only |
| `_prior_anomaly_fallback(p, precomp, reach)` | Fallback: prior medians + WSE anomaly scaling |
| `_wse_width_depth_correction(precomp, z0_post, n_post, priors)` | Post-hoc z₀ correction from WSE-Width regression (rectangular mode) |
| `initialize_theta(priors, precomp, months)` | Initial θ from prior medians |
| `laplace_uncertainty(obj, θ_map)` | Hessian-based posterior covariance |
| `safe_logpdf(d, x)` | ForwardDiff-safe log-pdf with smooth barriers |

### `preprocess.jl`
| Function | Description |
|---|---|
| `preprocess(xobs, Hobs, Wobs, Sobs; dx, min_slope, nknots, lambda)` | Build SWOTReach from raw observations |
| `drop_unobserved(x, H, W, S)` | Remove nodes with no valid observations |
| `build_chainage(xobs, dx)` | Construct uniform downstream-to-upstream chainage |
| `estimate_bed_slope(Sobs, min_slope)` | Time-averaged bed slope from SWOT slope data (uses abs(slope) to handle negative SWOT convention) |
| `obs_chainage(xobs, Hobs, Wobs, min_slope)` | Fill gaps in WSE/width via PCHIP, enforce monotonicity |
| `smooth_profile_bspline(x, y; nknots, lambda, monotone_slope)` | Penalized B-spline WSE smoothing |
| `gcv_select(x, y, B; lambda_range)` | GCV smoothing parameter selection |
| `slope_from_spline(spl, x)` | Slope from B-spline derivative |
| `enforce_monotonicity(spl, x, min_slope)` | Enforce minimum slope on B-spline profile |

### `priors.jl`
| Function | Description |
|---|---|
| `priors(sosfile, hmin, reachid; S0)` | Build priors from SoS database (monthly Q if available) |
| `priors(qwbm, hmin, class; reach)` | Build uninformative priors when SoS unavailable |
| `monthly_q_prior(priors, month)` | Get Q prior for a given calendar month |
| `z0_prior(qwbm, hmin, reach)` | Estimate z₀ prior from Manning depth scaling |
| `_manning_depth_estimate(q, S; n_est)` | Manning depth for z₀ prior estimation |
| `_build_q_prior(f, i, q_m, q_l, q_u)` | Build Q prior from SoS (monthly or single) |
| `_monthly_distributions(monthly_means, q_l, q_u, q_m)` | 12 monthly truncated LogNormals |
| `_single_q_prior(q_m, q_l, q_u)` | Single time-invariant truncated LogNormal |

---

## Inference Pipeline Flow

1. **Preprocess** (`preprocess`): Raw SWOT → `SWOTReach` (interpolation, bed slope, bankfull geometry)
2. **Priors** (`priors`): SoS database → `SWOTPriors` (monthly Q, n, r, z₀)
3. **Pre-compute** (`precompute_manning`): Extract geometry and observations for forward model; auto-estimate σ_obs and σ_W
4. **Initialize** (`initialize_theta`): Set θ₀ from prior medians (monthly Q, n=0.03, r=prior median, z₀=prior median)
5. **Optimize** (`infer` → L-BFGS on `neg_log_joint`):
   - Auto-select σ_obs from precomp if not specified
   - Adaptive λ_smooth: reduce for data-sparse reaches (data_fraction < 0.1)
   - If optimizer produces unphysical parameters OR Q spikes > 10× prior median: double σ_obs and retry
   - If after retry n > 0.05 (prior truncation) AND σ_obs was increased: fall back to prior+WSE anomaly
6. **Profile refit** (`_profile_refit`): Fix n and r from Step 5 at their MAP values, re-optimize only (Q, z₀). This breaks the Q–n degeneracy by constraining shape parameters, allowing the WSE observations to directly determine Q magnitude.
7. **WSE-Width depth correction** (`_wse_width_depth_correction`, rectangular only): After profile refit, regress observed W vs WSE per timestep to estimate dW/dH. If dW/dH > 0 with R² > 0.05, estimate depth range Δy = ΔW/dWdH and compute a corrected z₀. Only applies when z₀_post is significantly above the width-derived z₀ maximum. Often returns `nothing` for flat rivers where backwater dominates WSE.
8. **Uncertainty** (`laplace_uncertainty`): Invert Hessian at MAP → posterior covariance → Q_std via delta method

---

## Important Conventions

- **x coordinates**: downstream (x=0) to upstream (x increasing). This is the SWOT convention used throughout.
- **Downstream boundary condition**: `H_bc = reach.H[1, t]`. This is an INPUT to the forward model, NOT an observation. The likelihood only applies at nodes j ≥ 2.
- **Missing data**: SWOT observations use `missing` values extensively. Always check with `!ismissing()` or `!isnan()` (after conversion).
- **Depth positivity**: always clamp `y ≥ 0.01m` to prevent log/diverge issues.
- **Log-space for Q**: Q is always estimated in log-space (`logQ = θ[1:nt]`) to enforce Q > 0 and improve optimization conditioning.
- **Log-space for n and r**: Manning n and Dingman r are estimated as `logn = θ[nt+1]`, `logr = θ[nt+2]` with `n = exp(logn)`, `r = exp(logr)`.
- **z₀ in natural space**: z₀ = θ[nt+3] — not log-transformed, can be negative.
- **Yb depends on z₀**: `Yb_nodes = max.(hbf_z_nodes .- z₀, 0.01)` — the bankfull depth varies with the estimated bed elevation. The `hbf_z_nodes` field stores the z₀-independent part of Yb.

---

## Dependencies (Project.toml)

```toml
BSplineKit = "093aae92-e908-43d7-9660-e50ee39d5a0a"     # B-spline smoothing
DataInterpolations = "82cc6244-b520-54b8-b5a6-8a565e85f1d0" # PCHIP, LinearInterpolation
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"       # AD for gradient + Hessian
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
NCDatasets = "85f8d34a-cbdd-5861-8df4-14fed0d494ab"       # SoS/SWOT NetCDF I/O
Optim = "429524aa-4258-5aef-a3af-852621145aeb"             # L-BFGS optimizer
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
```

---

## Testing Strategy

### Unit Tests (`test/runtests.jl`)

1. **preprocess**: SWOTReach construction, chainage, bed slope, bankfull geometry, valid timesteps
2. **priors**: SWOTPriors structure, Q/n/r/z prior ranges, sampling
3. **cross-section geometry**: Dingman area, WSE, width calculations
4. **manning**: manning_Q/depth consistency, ForwardDiff compatibility, backwater length, boundary conditions, multi-parameter gradients
5. **infer**: precompute_manning, neg_log_joint (finite, gradient finite), basic inference convergence, NSE benchmarking, inference speed (<30s), safe_logpdf with ForwardDiff
6. **missing data**: all-H-missing, all-W-missing, drop_unobserved, degenerate reaches
7. **B-spline smoothing**: design matrix, penalty matrix, GCV, profile smoothing, slope extraction, monotonicity enforcement

### SWOT Reach Regression Tests (`test/swot_reachtests.jl`)

Located in `test/swot_reachtests.jl`. These tests use real SWOT data for three reaches to catch performance regressions. They are **skipped automatically** if data files are not found (set `SWOT_DATA` env var to the data root directory).

| Reach ID | Region | Notes |
|----------|--------|-------|
| 12291300071 | Africa (GRDC 1291100) | r ≥ 0.80, NSE ≥ 0.50, KGE ≥ 0.40 vs. SVS gage |
| 81270600061 | North America (USGS 4103820) | Q median ≤ 200 m³/s |
| 56395000161 | Oceania (small river) | σ_obs retry with wider bounds |

Run with:
```
cd Sad.jl && julia --project=. test/swot_reachtests.jl
```

---

## σ_obs Retry and Fallback Logic

The `infer` function includes a multi-stage convergence check:

1. **Initial optimization** with auto-estimated σ_obs (min floor 0.3 m)
2. **Physical parameter check**: n ∈ [0.005, 0.12], r ≤ 20, z₀ ∈ [hmin−30, hmin+10]
3. **Q spike check**: any Q > 10× monthly prior median
4. **n outside prior check**: n > 0.10 (unphysical)
5. **If any check fails**: double σ_obs and re-optimize (up to ~3 doublings)
6. **After retries**: if n > 0.05 AND σ_obs was increased → prior+WSE anomaly fallback
7. **Fallback** (`_prior_anomaly_fallback`): Q = Q_prior × exp(0.5 × 1.5 × log(depth_ratio)), clamped to [0.2×, 5×] prior, with uncertainty = max(prior_std, 0.3×Q)

### Adaptive λ_smooth

For data-sparse reaches (data_fraction < 0.1), λ_smooth is linearly scaled down to let monthly Q priors drive the seasonal cycle:
```
λ_use = λ_smooth × min(1, data_fraction / 0.1)
```

### Data fraction check

If less than 5% of timesteps are valid, inference is skipped entirely and the result is flagged as `fallback=true` with Q = prior median + NaN uncertainty.

---

## Resolved Issues

### σ_obs Retry for Small Rivers
Reach 56395000161 (Oceania, mean Q ≈ 19 m³/s) has extremely precise WSE observations (median noise ≈ 0.0 m), auto-estimated σ_obs ≈ 0.15 m (floored to 0.3 m). The Manning+Dingman model has structural error ≫ 0.3 m. The retry mechanism automatically finds σ_obs = 1.2 m (4× initial), producing n=0.071, r=5.48.

### n Prior Truncation at 0.05
The n prior is truncated to [0.01, 0.05] as regularisation. Values above 0.05 typically indicate the model is compensating for structural error. When σ_obs retry is triggered AND the final n exceeds 0.05, inference falls back to the prior+WSE anomaly method, which avoids overfitting the WSE data with an excessively rough n.

### Q Spike Detection
Any timestep where Q > 10× the monthly prior median triggers the σ_obs retry. This catches degenerate solutions where the optimizer inflates Q to compensate for poor channel parameter estimates.

### SWOT Slope Sign Convention
SWOT water surface slope is reported in oceanographic convention (negative for normally-flowing rivers where WSE decreases downstream). The `estimate_bed_slope` function uses `abs(slope)` to handle this correctly. About 6% of reaches have ALL negative SWOT slopes; previously these fell to the `min_slope = 1e-5` floor, producing S₀ 12–85× too small. This caused the optimizer to compensate with unphysical Manning n > 0.05.

### Q-z₀ Degeneracy in Rectangular Mode
In rectangular channel mode (W = Wb constant), WSE = z₀ + y but width provides no independent depth constraint. The optimizer pushes z₀ toward hmin (shallow depth, low Q). Mitigations: (1) truncated Normal z₀ prior centred at Manning depth estimate with upper bound < hmin, (2) stiffer barrier penalty at truncation bounds, (3) WSE-Width depth correction using minimum observed width as bankfull reference. For flat rivers (S₀ < 1e-3) with strong backwater, the downstream boundary condition controls more WSE variance than local Q, creating a fundamental limitation of the 1D steady-state Manning model.