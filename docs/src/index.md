# SAD.jl Documentation

```@contents
Pages = ["installation.md", "algorithm.md", "use_cases.md", "api.md"]
```

The Surface Water Ocean Topography (SWOT) satellite mission, launched in December 2022, provides unprecedented observations of river water surface elevation (WSE), width, and slope at global scale. Because river discharge is not directly observed by SWOT, the SAD (SWOT Assimilated Discharge) algorithm estimates discharge from these observables using a data assimilation framework.

**SAD v2** implements a **joint MAP (Maximum A Posteriori) estimation** approach that simultaneously infers discharge, Manning roughness, Dingman channel shape, and bed elevation by combining an **analytical forward model** with prior distributions from the SWORD of Science (SoS) database.

Key features of v2:

- **Analytical forward model**: Manning equation with first-order GVF (Gradually Varied Flow) backwater perturbation — no ODE solve, fully AD-compatible
- **Joint optimization**: All parameters (Q₁…Qₜ, n, r, z₀) estimated simultaneously via L-BFGS, eliminating the two-stage bias of v1
- **Monthly Q priors**: Seasonal discharge priors from SoS, resolving n–z₀–Q equifinality
- **Robust likelihood**: Student-t (ν ≈ 3–5) for outlier resistance
- **Temporal smoothness**: Penalty on successive log-Q differences for coherent hydrographs
- **Laplace uncertainty**: Posterior covariance from the Hessian at MAP
- **Adaptive σ_obs**: Auto-estimated WSE noise with retry mechanism for small rivers
- **Fallback strategy**: Prior + WSE anomaly when MAP produces unphysical parameters