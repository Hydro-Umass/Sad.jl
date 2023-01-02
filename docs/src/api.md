# Library

## Types

```@docs
	Sad.Dingman
	Sad.Rectangular
```

## Functions

### Cross sections

```@docs
	Sad.width
	Sad.depth
```

### Priors

```@docs
	Sad.get_samples
	Sad.lhs_ensemble
	Sad.prior_ensemble
	Sad.gvf_ensemble!
	Sad.rejection_sampling
	Sad.priors
```

### Kalman filters

```@docs
	Sad.letkf
```

### Gradually Varied Flow

```@docs
	Sad.froude
	Sad.interpolate_hbc
	Sad.dydx
	Sad.gvf
```

### Sad

```@docs
	Sad.drop_unobserved
	Sad.assimilate
	Sad.bathymetry!
	Sad.ahg_constrain!
	Sad.calc_bed_slope
	Sad.calc_dA
	Sad.flow_parameters
	Sad.estimate
```
