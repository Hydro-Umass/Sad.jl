# Use cases

!!! warning "Array ordering in SAD"
    SAD expects time series of hydraulic variable profiles as 2-D arrays, with
    the **row dimension** being the cross-sections from **downstream (index 1) to
    upstream (last index)** and the column dimension being time. Missing data are
    acceptable — use [Julia's `missing` type](https://docs.julialang.org/en/v1/manual/missing/)
    rather than `NaN`.

## Basic usage with synthetic data (Pepsi)

Load a river reach from a NetCDF file and preprocess the SWOT observations:

```julia
using NCDatasets, Statistics, Sad

f = Dataset("data/pepsi1/Po.nc")
g = NCDatasets.group(f, "XS_Timeseries")
qwbm = Float64(NCDatasets.group(f, "River_Info")["QWBM"][1])
x = (g["X"][:][end] .- g["X"][:])[end:-1:1, 1]
H = convert(Matrix{Sad.FloatM}, g["H"][end:-1:1, :])
W = convert(Matrix{Sad.FloatM}, g["W"][end:-1:1, :])

# Compute slope from WSE
S = diff(H, dims=1) ./ diff(x)
S = convert(Matrix{Sad.FloatM}, vcat(S[1:1, :], S))

# Preprocess observations into a SWOTReach
reach = Sad.preprocess(x, H, W, S)
```

Construct priors. When SoS data is available, monthly Q priors are automatically
constructed:

```julia
# From SoS database (recommended)
sosfile = "data/sos/na_sword_v17_SOS_priors.nc"
p = Sad.priors(sosfile, reach.hmin, reachid; S0=mean(reach.S0.(reach.x)))
```

When SoS data is unavailable, construct uninformative priors from the basin
mean discharge and river planform class:

```julia
# Uninformative priors
p = Sad.priors(qwbm, reach.hmin, Sad.sinuous; reach=reach)
```

Run SAD v2 inference:

```julia
# With monthly Q priors from observation timestamps
time_str = ["2023-01-15", "2023-02-12", ...]  # one per timestep
result = Sad.infer(p, reach; time_str=time_str)
```

Or with explicit keyword arguments for advanced control:

```julia
result = Sad.infer(p, reach;
    time_str    = time_str,     # timestamps for monthly Q-prior selection
    σ_obs       = NaN,          # NaN = auto-estimate from data
    ν           = 5.0,         # Student-t degrees of freedom (Inf = Gaussian)
    λ_smooth    = 0.1,          # temporal smoothness weight on log-Q
    iterations  = 500,          # max L-BFGS iterations
    g_tol       = 1e-6,         # gradient convergence tolerance
    use_width   = false,        # include width likelihood
)
```

The result is a `NamedTuple` with:

| Field | Type | Description |
|-------|------|-------------|
| `Q_post` | `Vector{Float64}` | Posterior mean Q per timestep (MAP) |
| `Q_std` | `Vector{Float64}` | Posterior std of Q (Laplace delta method) |
| `n_post` | `Float64` | MAP Manning n |
| `r_post` | `Float64` | MAP Dingman r |
| `z0_post` | `Float64` | MAP downstream bed elevation |
| `θ_map` | `Vector{Float64}` | Full MAP parameter vector |
| `Σ` | `Matrix{Float64}` | Posterior covariance (Laplace) |
| `result` | `Optim.OptimResult` | Optimization result |
| `valid_ts` | `Vector{Int}` | Indices of valid timesteps |
| `precomp` | `ManningPrecomp` | Pre-computed data |
| `fallback` | `Bool` | Whether prior+WSE anomaly fallback was used |

### Extracting and plotting results

```julia
# Discharge time series
Q = result.Q_post
Q_lo = result.Q_post .- 1.96 .* result.Q_std
Q_hi = result.Q_post .+ 1.96 .* result.Q_std

# Physical parameters
println("Manning n = $(round(result.n_post, digits=4))")
println("Dingman r = $(round(result.r_post, digits=2))")
println("Bed elevation z₀ = $(round(result.z0_post, digits=2)) m")

# Check if fallback was used
if result.fallback
    @warn "Prior + WSE anomaly fallback was used (n likely exceeded prior truncation)"
end
```

## Confluence: operational SWOT processing

When running SAD within the Confluence framework, priors are obtained from the
SWORD network database and SoS priors file:

```julia
swordfile = "data/sword/na_sword_v17.nc"
sosfile   = "data/sos/na_sword_v17_SOS_priors.nc"
swotfile  = "data/swot/74267400111_SWOT.nc"

# Read SWOT observations
nids, x = Sad.river_info(reachid, swordfile)
H, W, S, dA, Hr, Wr, Sr, time_str = Sad.read_swot_obs(swotfile, nids)
H, W, S = Sad.drop_unobserved(H, W, S)  # remove fully-missing nodes

# Preprocess and build priors
reach = Sad.preprocess(x, H, W, S)
p = Sad.priors(sosfile, reach.hmin, reachid; S0=mean(reach.S0.(reach.x)))

# Run inference with monthly priors
result = Sad.infer(p, reach; time_str=time_str)
```

## Handling small rivers and σ_obs failures

For small rivers where the auto-estimated σ_obs is much smaller than the model
structural error, the σ_obs retry mechanism automatically doubles σ_obs until
convergence. No user intervention is needed — the retry logic is built into
`infer`.

However, you may want to override σ_obs manually for debugging:

```julia
# Force a specific σ_obs (meters)
result = Sad.infer(p, reach; σ_obs=2.0)

# Disable temporal smoothness
result = Sad.infer(p, reach; λ_smooth=0.0)

# Use Gaussian likelihood instead of Student-t
result = Sad.infer(p, reach; ν=Inf)
```

## Missing data handling

SAD v2 handles missing SWOT observations gracefully. In the observation matrix,
use `missing` values for missing data:

```julia
using Sad
H = convert(Matrix{Sad.FloatM}, H_raw)  # FloatM = Union{Missing, Float64}
W = convert(Matrix{Sad.FloatM}, W_raw)
```

During preprocessing, nodes with no valid observations across all timesteps are
removed via `drop_unobserved`. During inference, the likelihood automatically
skips nodes with missing WSE or width observations. Timesteps with fewer than
2 valid WSE observations are flagged as invalid and excluded from the MAP
objective.