# SWOT Assimilated Discharge (SAD)

**Documentation** [![](https://img.shields.io/badge/docs-online-blue.svg)](https://hydro-umass.github.io/Sad.jl/)

**Build/Tests** [![CI](https://github.com/Hydro-Umass/Sad.jl/workflows/CI/badge.svg?branch=master)](https://github.com/Hydro-Umass/Sad.jl/actions?query=workflow:CI)

The Surface Water Ocean Topography (SWOT) satellite mission that launched in December 2022 will offer a unique opportunity to map river discharge at an unprecedented spatial resolution globally from observations of water surface elevation, width, and slope. Because river discharge will not be directly observed from SWOT, a number of algorithms of varying complexity have been developed to estimate discharge from SWOT observables. The SAD algorithm operates on the set of SWOT observables (i.e., WSE, width, and slope) and derives an estimate of river discharge and its associated uncertainty using a data assimilation scheme. The assimilation scheme involves the “first‐guess” estimation of hydraulic variables by combining a forward model with a set of prior probability distributions before assimilating the SWOT observations. The priors are estimated with a sampling approach, as these data‐driven methods have shown promise in many fields. The objective of the algorithm is to estimate discharge at each river reach when SWOT observations become available.
