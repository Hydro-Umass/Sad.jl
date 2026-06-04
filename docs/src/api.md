# Library

## Types

### SWOTObs

```@docs
Sad.SWOTObs
```

### SWOTReach

```@docs
Sad.SWOTReach
```

### SWOTPriors

```@docs
Sad.SWOTPriors
```

### ManningPrecomp

```@docs
Sad.ManningPrecomp
```

### River

```@docs
Sad.River
```

## Manning forward model

These functions implement the analytical Manning + first-order GVF perturbation
forward model. They are fully AD-compatible (work with `ForwardDiff.Dual` numbers).

```@docs
Sad.manning_Q
Sad.manning_depth
Sad.backwater_length
Sad.backwater_depth
Sad.manning_wse
Sad.manning_wse_backwater
```

### Dingman cross-section geometry

Utility functions for the Dingman power-law cross-section.

```@docs
Sad.area
Sad.compute_wse
Sad.compute_width
```

## Preprocessing

### Main entry point

```@docs
Sad.preprocess
```

### Pre-computation for inference

```@docs
Sad.precompute_manning
```

### Helper functions

```@docs
Sad.build_chainage
Sad.estimate_bed_slope
Sad.drop_unobserved
Sad.obs_chainage
Sad.interpolate_to_chainage
```

### B-spline smoothing

Penalized cubic B-spline functions for WSE profile smoothing and slope
estimation. These replace PCHIP for robust reach-averaged slope extraction.

```@docs
Sad.smooth_profile_bspline
Sad.gcv_select
Sad.spline_design_matrix
Sad.penalty_matrix
Sad.slope_from_spline
Sad.enforce_monotonicity
```

## Priors

### From SoS database (recommended)

Construct priors from the SWORD of Science database, including monthly
climatological discharge priors when available.

```@docs
Sad.priors
```

See also: [`Sad.monthly_q_prior`](@ref) for monthly Q prior access.

### Monthly Q prior

```@docs
Sad.monthly_q_prior
```

### Bed elevation prior

```@docs
Sad.z0_prior
```

## Inference

### Main entry point

```@docs
Sad.infer
```

The `infer` function runs the complete SAD v2 pipeline:

1. Pre-compute reach geometry and observations (`precompute_manning`)
2. Initialize parameters from prior medians (`initialize_theta`)
3. L-BFGS optimization of `neg_log_joint` with ForwardDiff gradients
4. σ_obs retry for unphysical parameters
5. Prior + WSE anomaly fallback if n > 0.05 after retry
6. Laplace uncertainty from the Hessian at MAP

### Objective function

```@docs
Sad.neg_log_joint
```

### Uncertainty

```@docs
Sad.laplace_uncertainty
```

### Noise estimation

```@docs
Sad.estimate_σ_obs
Sad.estimate_σ_W
```

### Safe log-pdf

```@docs
Sad.safe_logpdf
```

## I/O helpers

The `read_swot_obs` and `river_info` functions are defined in
`scripts/swot.jl` (not exported from the module). They are used for loading
SWOT observation data from NetCDF files in operational processing scripts.

## Type alias

```@docs
Sad.FloatM
```