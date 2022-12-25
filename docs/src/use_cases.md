# Use cases

## Pepsi-1 experiment

We start by loading the data, including river node information and synthetic SWOT observations

```julia
    using NCDatasets, Statistics, Distributions, Sad
    f = Dataset("../../data/pepsi1/Po.nc")
    g = NCDatasets.group(f, "XS_Timeseries")
	qwbm = NCDatasets.group(f, "River_Info")["QWBM"][1]
    x = (g["X"][:][end] .- g["X"][:])[end:-1:1, 1]
    Q = g["Q"][:][end:-1:1, :]
    H = g["H"][:][end:-1:1, :]
    W = g["W"][:][end:-1:1, :]
```

Bankfull width and water surface elevation can be guessed as the maximum from the observed time series, while an initial estimate of bed slope can be obtained from the mean of surface water slope

```julia
    hbf = maximum(H, dims=2)[:, 1]
    wbf = maximum(W, dims=2)[:, 1]
    S0 = mean(diff(H, dims=1) ./ diff(x), dims=2)[:, 1]
    S0 = [S0[1]; S0]
```

Then we can derive the prior distributions using rejection sampling from uninformative priors

```julia
	Qp0, np0, rp0, zp0 = Sad.priors(qwbm, minimum(H[1, :]))
	Qp, np, rp, zp = Sad.rejection_sampling(Qp0, np0, rp0, zp0, x, H, S0, mean(H[1, :]), wbf, hbf, 1000);
```

The ensemble of discharge, roughness coefficient, channel shape parameter and bed elevation can now be generated

```julia
	Qe, ne, re, ze = Sad.prior_ensemble(x, Qp, np, rp, zp, 1000);
```

Bed elevation and slope are then estimated by assimilating the time-average water surface elevation profile

```julia
	S = diff(H, dims=1) ./ diff(x);
    Se = repeat(S', outer=3)'[:, 1:1000];
    Se = [Se[1, :]'; Se]
    Sad.bathymetry!(ze, Se, Qe, ne, re, x, hbf, wbf, mean(H, dims=2)[:, 1])
	zp = Truncated(Normal(mean(ze[1, :]), std(ze[1, :])), -Inf, minimum(H[1, :]))
```

and finally we assimilate the observed water surface elevation to estimate discharge and flow parameters, which include roughness coefficient and minimum cross-sectional area

```julia
	Qa, na = assimilate(H, W, x, wbf, hbf, Se, Qp, np, rp, zp, nens)
```

## Confluence
