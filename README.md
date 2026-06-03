# SWOT Assimilated Discharge (SAD) v2

**Documentation** [![](https://img.shields.io/badge/docs-online-blue.svg)](https://hydro-umass.github.io/Sad.jl)

**Build/Tests** [![CI](https://github.com/Hydro-Umass/Sad.jl/workflows/CI/badge.svg)](https://github.com/Hydro-Umass/Sad.jl/actions?query=workflow:CI)

**Code coverage**
[![codecov](https://codecov.io/github/Hydro-Umass/Sad.jl/branch/main/graph/badge.svg?token=80Z2DDOA8O)](https://codecov.io/github/Hydro-Umass/Sad.jl)

The Surface Water Ocean Topography (SWOT) satellite mission provides observations of water surface elevation (WSE), width, and slope at an unprecedented spatial resolution globally. Because river discharge is not directly observed, the SAD algorithm estimates discharge by combining a Manning-based hydraulic forward model with prior distributions in a maximum a posteriori (MAP) inference framework. SAD v2 replaces the LETKF assimilation scheme from v1 with deterministic optimisation (L-BFGS), offering faster runtimes and improved handling of sparse or noisy SWOT observations.

## Overview

SAD v2 estimates discharge for individual river reaches by:

1. **Preprocessing** SWOT observations of WSE, width, and slope, and deriving reach geometry (bankfull width, bed elevation, slope) from the SWOT timeseries and [SWORD](https://github.com/SWOT-Data/SWORD) network data.
2. **Computing priors** for Manning's roughness coefficient (*n*), channel shape parameter (*r*), bed elevation (*z₀*), and monthly discharge from the [SoS](https://doi.org/10.1029/2023WR036123) database or defaults.
3. **MAP inference** — maximising the joint log-posterior of log-discharge (per timestep), log-roughness, log-shape, and bed elevation using L-BFGS, with the Manning + Dingman forward model as the likelihood and informative priors as regularisation.

### Manning + Dingman forward model

Water depth at each cross-section is computed via the gradually varied flow (GVF) equation:

$$\frac{dY}{dx} = \frac{S_0 - S_f}{1 - \mathrm{Fr}^2}$$

where friction slope *S_f* follows Manning's equation and the Froude number uses the Dingman cross-section geometry with shape parameter *r*. The model predicts WSE and width at observed nodes, which are compared against SWOT observations in the likelihood.

### Priors

| Parameter | Prior | Source |
|-----------|-------|--------|
| Manning's *n* | Truncated Normal(0.03, 0.01; [0.01, 0.05]) | SoS database, regularised |
| Shape *r* | Truncated Normal from SoS | SoS database |
| Bed elevation *z₀* | Uniform(hmin − depth, hmin + depth) | Manning depth estimate |
| Discharge *Q* | Truncated LogNormal (monthly) | SoS monthly ML priors |

### Uncertainty estimation

SAD v2 provides two methods for discharge uncertainty, depending on the inference path:

**Laplace approximation (MAP inference).** When the optimiser converges to a physical MAP estimate, posterior uncertainty is computed via the Laplace (curvature) approximation. The Hessian of the negative log-posterior is evaluated at the MAP point using automatic differentiation (ForwardDiff), then inverted to obtain the approximate posterior covariance matrix Σ. Discharge uncertainty is derived via the delta method on the log-Q parameterisation:

$$\sigma_Q \approx Q \cdot \sigma_{\log Q}$$

where σ_{log Q} = √(Σ_{ii}) for the i-th timestep. This captures both observation uncertainty (through the likelihood) and parameter uncertainty (through the prior curvature and parameter correlations). Negative variances from non-convexity are handled by adding a small regularisation term to the Hessian diagonal.

**Prior-based uncertainty (WSE anomaly fallback).** When MAP inference fails and the algorithm falls back to the WSE anomaly method, the prior standard deviation is used as the uncertainty baseline. For each month, the prior IQR (interquartile range / 1.35) provides a robust standard deviation estimate. If the monthly prior is too wide or degenerate, a fallback using the prior median as the standard deviation is used. The reported uncertainty is the maximum of the prior standard deviation and 30% of the estimated Q, ensuring that fallback uncertainties are always at least as large as the prior uncertainty.

Both methods produce per-timestep discharge standard deviations (`Q_std` in the return value and `Q_u` in the output NetCDF).

### Robustness features

SAD v2 includes several guards for real-world SWOT data:

- **Data-sparse reaches** (< 5% valid timesteps) return invalid/NaN rather than unreliable estimates.
- **Adaptive σ_obs** — WSE observation noise is auto-estimated from cross-section residuals and scaled by data fraction.
- **Convergence retry** — when the optimiser converges to unphysical parameters (e.g., *n* > 0.12) or produces Q spikes (> 10× the monthly prior median), σ_obs is progressively doubled until physical parameters are obtained.
- **WSE anomaly fallback** — when σ_obs retry was required *and* the final Manning's *n* exceeds the prior truncation (0.05), the algorithm falls back to a depth-scaling anomaly method: *Q = Q_prior × exp(α · log(Y/Y₀))*, which avoids bias from an excessively rough *n* compensating for model structural error.
- **NaN/bounds guards** in priors protect against degenerate Uniform(a, b) distributions when hmin = 0 or slope is missing.

## Quick start

```julia
using Sad, NCDatasets, Statistics, Distributions

include("scripts/swot.jl")

# Read reach data
reachid = 12291300071
nids, x = river_info(reachid, "data/sword/af_sword_v17.nc")
H, W, S, dA, Hr, Wr, Sr, time_str = read_swot_obs("data/swot/$(reachid)_SWOT.nc", nids)

# Preprocess and compute priors
reach = Sad.preprocess(x, H, W, S)
p = Sad.priors("data/sos/af_sword_v17_SOS_priors.nc", reach.hmin, reachid;
                S0=mean(reach.S0.(reach.x)))

# Run MAP inference
res = Sad.infer(p, reach; time_str=time_str, λ_smooth=0.1)

# Results
println("n = $(res.n_post), r = $(res.r_post), z0 = $(res.z0_post)")
println("Q range: $(minimum(res.Q_post)) – $(maximum(res.Q_post))")
```

## Running the full pipeline

The `scripts/swot.jl` script provides `main(reachid)` and `write_output()` for batch processing of SWOT reaches. The `run_all.jl` script processes a set of science test reaches.

```bash
julia --project=. run_all.jl
```

Output NetCDF files are written with discharge estimates (`Qa`, `Q_u`), Manning's *n*, shape *r*, bed elevation *z₀*, and metadata including `valid`, `fallback`, and `converged` flags.

## Installation

```julia
import Pkg
Pkg.add("Sad")
Pkg.instantiate()
```

## Running tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## References

- Andreadis, K. M., Brinkerhoff, C. B., & Gleason, C. J. (2020). Constraining the assimilation of SWOT observations with hydraulic geometry relations. *Water Resources Research*, 56, e2019WR026811. https://doi.org/10.1029/2019WR026811

- Dingman, S. L., & Afshari, S. (2018). Field verification of analytical at-a-station hydraulic-geometry relations. *Journal of Hydrology*, 564, 859–872. https://doi.org/10.1016/j.jhydrol.2018.07.020
