# Installation instructions

The latest version of Sad.jl can be installed with Julia's [built-in package manager](https://docs.julialang.org/en/v1/stdlib/Pkg/). You can add the package and then run `instantiate` to install dependencies

```julia
import Pkg
Pkg.add("Sad")
Pkg.instantiate()
```

Similarly the package can be updated 

```julia
Pkg.update("Sad")
```

!!! warn "Use Julia 1.6 or newer"
	The Sad.jl code has only been tested with Julia v1.6 and above, so please update your installation accordingly.
