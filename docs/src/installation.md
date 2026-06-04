# Installation

The latest version of SAD.jl can be installed with Julia's built-in [package manager](https://docs.julialang.org/en/v1/stdlib/Pkg/):

```julia
import Pkg
Pkg.add("Sad")
Pkg.instantiate()
```

Or from the repository source:

```julia
import Pkg
Pkg.add(url="https://github.com/Hydro-Umass/Sad.jl")
Pkg.instantiate()
```

To update:

```julia
Pkg.update("Sad")
```

## Julia version

!!! warn "Julia 1.10 or newer required"
    SAD.jl v2 requires Julia 1.10 or later due to dependencies on ForwardDiff and Optim.

## Dependencies

SAD.jl v2 depends on the following packages:

| Package | Purpose |
|---------|---------|
| [Optim.jl](https://github.com/JuliaNlsolvers/Optim.jl) | L-BFGS optimizer for MAP estimation |
| [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) | Automatic differentiation for gradients and Hessians |
| [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) | Prior distributions (LogNormal, truncated Normal, etc.) |
| [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl) | PCHIP and linear interpolation for SWOT profiles |
| [BSplineKit.jl](https://github.com/jw3126/BSplineKit.jl) | B-spline smoothing for WSE profiles |
| [NCDatasets.jl](https://github.com/Alexander-Barth/NCDatasets.jl) | Reading SWOT and SoS NetCDF data |

The v2 release removes the Dependency on DifferentialEquations.jl (previously used for the GVF ODE solver), as the forward model is now fully analytical.

## Testing

Run the test suite with:

```bash
cd Sad.jl && julia --project=. -e 'using Pkg; Pkg.test()'
```

### SWOT reach regression tests

Regression tests using real SWOT data require large data files. Set the `SWOT_DATA` environment variable to point to the data directory:

```bash
cd Sad.jl && julia --project=. test/swot_reachtests.jl
```