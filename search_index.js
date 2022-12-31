var documenterSearchIndex = {"docs":
[{"location":"api/#Library","page":"API","title":"Library","text":"","category":"section"},{"location":"api/#Types","page":"API","title":"Types","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"\tSad.Dingman\n\tSad.Rectangular","category":"page"},{"location":"api/#Sad.Dingman","page":"API","title":"Sad.Dingman","text":"Dingman river channel cross section.\n\nUses generalized expressions for cross-section geometry based on at-a-station hydraulic geometry relations.\n\nWb: bankfull width\nYb: bankfull depth\nYm: flow depth\nr: geometry exponent\nS0: bed slope\nn: Manning roughness coefficient\n\nSee also: Rectangular\n\n\n\n\n\n","category":"type"},{"location":"api/#Sad.Rectangular","page":"API","title":"Sad.Rectangular","text":"Rectangular river channel cross section.\n\nWb: bankfull width\nYb: bankfull depth\nYm: flow depth\nS0: bed slope\nn: Manning roughness coefficient\n\nSee also: Dingman\n\n\n\n\n\n","category":"type"},{"location":"api/#Functions","page":"API","title":"Functions","text":"","category":"section"},{"location":"api/#Cross-sections","page":"API","title":"Cross sections","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"\tSad.width\n\tSad.depth","category":"page"},{"location":"api/#Sad.width","page":"API","title":"Sad.width","text":"width(xs)\n\nCalculate water surface width of cross section.\n\nArguments\n\nxs: cross section\n\n\n\n\n\n","category":"function"},{"location":"api/#Sad.depth","page":"API","title":"Sad.depth","text":"depth(xs)\n\nCalculate average depth of cross section.\n\nArguments\n\nxs: cross section\n\n\n\n\n\n","category":"function"},{"location":"api/#Priors","page":"API","title":"Priors","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"\tSad.get_samples\n\tSad.lhs_ensemble\n\tSad.prior_ensemble\n\tSad.gvf_ensemble!\n\tSad.rejection_sampling\n\tSad.priors","category":"page"},{"location":"api/#Sad.get_samples","page":"API","title":"Sad.get_samples","text":"get_samples(p, samples)\n\nGet quantiles from p according to weights provided by Latin Hypercube Sampling samples.\n\n\n\n\n\n","category":"function"},{"location":"api/#Sad.lhs_ensemble","page":"API","title":"Sad.lhs_ensemble","text":"lhs_ensemble(nens, args...)\n\nGenerate an ensemble of size nens using Latin Hypercube Sampling of the list of distributions or collections provided as args.\n\n\n\n\n\n","category":"function"},{"location":"api/#Sad.prior_ensemble","page":"API","title":"Sad.prior_ensemble","text":"prior_ensemble(x, Qp, np, rp, zp, nens)\n\nGenerate a prior ensemble of discharge, roughness coefficient, channel shape parameter, and downstream bed elevation from provided distributions or sample data.\n\nArguments\n\nx: distances between cross sections\nQp: discharge prior distribution\nnp: roughness coefficient distribution\nrp: channel shape parameter distribution\nzp: downstream bathymetry distribution\nnens: ensemble size\n\n\n\n\n\n","category":"function"},{"location":"api/#Sad.gvf_ensemble!","page":"API","title":"Sad.gvf_ensemble!","text":"gvf_ensemble!(H::Vector{Float64}, S, x::Vector{Float64}, hbf::Vector{Float64}, wbf::Vector{Float64}, Qe::Vector{Float64}, ne::Vector{Float64}, re::Vector{Float64}, ze)\n\nGenerate an ensemble of water height profiles from Gradually-Varied-Flow simulations and associated profiles of bed elevation.\n\nArguments\n\nH: water surface elevation\nS: bed slope\nx: channel chainage\nhbf: bankfull water surface elevation\nwbf: bankfull width\nQe: ensemble discharge\nne: ensemble roughness coefficient\nre: ensemble channel shape parameter\nze: ensemble bed elevation profiles\n\n\n\n\n\n","category":"function"},{"location":"api/#Sad.rejection_sampling","page":"API","title":"Sad.rejection_sampling","text":"rejection_sampling(Qp, np, rp, zp, x, H, S0, hbc, wbf, hbf, nens, nsamples)\n\nUse rejection sampling to select a subset of the prior ensemble.\n\nArguments\n\nQp: prior discharge distribution\nnp: prior roughness coefficient distribution\nrp: prior channel shape parameter distribution\nzp: prior distribution for downstream bed elevation\nx: channel chainage\nH: time series of observed water surface elevation profiles\nS0: prior estimate of channel bed slope\nhbc: mean downstream flow depth used as boundary condition\nwbf: bankfull width\nhbf: bankfull water surface elevation\nnens: ensemble size\nnsamples: number of samples\n\n\n\n\n\n","category":"function"},{"location":"api/#Sad.priors","page":"API","title":"Sad.priors","text":"priors(ncfile)\n\nDerive prior distributions for discharge, roughness coefficient, channel shape parameter, and bed elevation by reading the SWORD a-priori database.\n\nArguments\n\nncfile: NetCDF file of SWORD database\n\n\n\n\n\npriors(qwbm, hmin)\n\nDerive distributions for discharge, roughness coefficient, channel shape parameter, and bed elevation using uninformative priors.\n\nArguments\n\nqwbm: mean discharge\nhmin: minimum downstream water surface elevation\n\n\n\n\n\n","category":"function"},{"location":"api/#Kalman-filters","page":"API","title":"Kalman filters","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"\tSad.letkf","category":"page"},{"location":"api/#Sad.letkf","page":"API","title":"Sad.letkf","text":"letkf(A, d, HA, E, ρ, diagR)\n\nApply the Local Ensemble Transform Kalman Filter algorithm.\n\nArguments\n\nA: state variable ensemble matrix\nd: observation vector\nHA: model-predicted observation ensemble matrix\nE: observation error ensemble matrix\nρ: covariance inflation parameter (optional, default value of 1.05)\ndiagR: force observation error covariance to be diagonal, i.e. no spatial correlation (optional, default is false)\n\n\n\n\n\n","category":"function"},{"location":"api/#Gradually-Varied-Flow","page":"API","title":"Gradually Varied Flow","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"\tSad.froude\n\tSad.interpolate_hbc\n\tSad.dydx\n\tSad.gvf","category":"page"},{"location":"api/#Sad.froude","page":"API","title":"Sad.froude","text":"froude(Q::Float64, xs::CrossSection)\n\nCalculate the Froude number.\n\nArguments\n\nQ: discharge\nxs: cross section\n\n\n\n\n\n","category":"function"},{"location":"api/#Sad.interpolate_hbc","page":"API","title":"Sad.interpolate_hbc","text":"interpolate_h(x, H, S)\n\nInterpolate water surface elevation boundary condition by assuming uniform flow, i.e., energy slope is equal to bed slope.\n\nArguments\n\nx: channel chainage\nH: water surface elevation profile\nS: bed slope\n\n\n\n\n\n","category":"function"},{"location":"api/#Sad.dydx","page":"API","title":"Sad.dydx","text":"dydx(y, p, x)\n\nOrdinary differential equation describing the Gradually-Varied-Flow model.\n\n\n\n\n\n","category":"function"},{"location":"api/#Sad.gvf","page":"API","title":"Sad.gvf","text":"Calculate a water surface profile by solving the Gradually-Varied-Flow equation.\n\nQ: discharge\nybc: downstream boundary condition for depth\nS0: bed slope for reach\nn: roughness coefficient\nx: downstream distance for each cross section\nwbf: bankfull width for each cross section\nybf: bankfull depth for each cross section\nr: channel geometry coefficient for each cross section\n\n\n\n\n\n","category":"function"},{"location":"api/#Sad","page":"API","title":"Sad","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"\tSad.drop_unobserved\n\tSad.assimilate\n\tSad.bathymetry!\n\tSad.ahg_constrain!\n\tSad.flow_parameters\n\tSad.estimate","category":"page"},{"location":"api/#Sad.drop_unobserved","page":"API","title":"Sad.drop_unobserved","text":"drop_unobserved(x, H, W)\n\nRemove cross sections with no valid observations.\n\nArguments\n\nx: channel chainage\nH: time series of water surface elevation profiles\nW: time series of width profiles\n\n\n\n\n\n","category":"function"},{"location":"api/#Sad.assimilate","page":"API","title":"Sad.assimilate","text":"assimilate(H, W, S, x, wbf, hbf, S0, Qp, np, rp, zp, nens, constrain)\n\nAssimilate SWOT observations for river reach.\n\nArguments\n\nH: time series of water surface elevation profiles\nW: time series of water surface width profiles\nS: time series of water surface slope profiles\nx: downstream distance for each cross section\nwbf: bankfull width\nhbf: bankfull depth\nQₚ: prior probability distribution for discharge\nnₚ: prior probability distribution for roughness coefficient\nrₚ: prior probability distribution for channel shape parameter\nzₚ: prior distribution for downstream bed elevation\nnens: ensemble size\nconstrain: switch for applying AHG constraint\n\n\n\n\n\n","category":"function"},{"location":"api/#Sad.bathymetry!","page":"API","title":"Sad.bathymetry!","text":"bathymetry!(ze, Se, Qe, ne, re, x, hbf, wbf, H, ϵₒ)\n\nEstimate channel bed slope and thalweg (elevation) by assimilating SWOT water surface elevation observations.\n\nArguments\n\nze: ensemble of bed elevations\nSe: ensemble of bed slopes\nQe: ensemble of discharge\nne: ensemble of roughness coefficient\nre: ensemble of channel shape parameters\nx: channel chainage\nhbf: bankfull water surface elevation\nwbf: bankfull width\nH: time-averaged water surface elevation observed\nϵₒ: observation error (default value of 10 cm)\n\n\n\n\n\n","category":"function"},{"location":"api/#Sad.ahg_constrain!","page":"API","title":"Sad.ahg_constrain!","text":"ahg_constrain!(Qa, hbc, hbf, wbf, S0, x, ne, re, ze, tol=1e-3)\n\nConstrain assimilated discharge with At-a-Station hydraulic geometry (AHG) relationships.\n\nArguments\n\nQa: assimilated ensemble discharge\nhbc: downstream boundary water surface elevation\nhbf: bankfull water surface elevation\nwbf: bankfull width\nS0: bed slope\nx: channel chainage\nne: ensemble roughness coefficient\nre: ensemble channel shape parameter\nze: ensemble bed elevation\ntol: numeric tolerance for closing the AHG relationships\n\n\n\n\n\n","category":"function"},{"location":"api/#Sad.flow_parameters","page":"API","title":"Sad.flow_parameters","text":"flow_parameters(Qa, na, x, H, W, S, S0, hbf, wbf, r, z)\n\nEstimate flow parameters (roughness coefficient and baseflow cross-sectional area) from assimilated discharge.\n\nArguments\n\nQa: time series of estimated discharge\nna: time series of estimate roughness coefficient\nx: channel chainage\nH: time series of observed water surface elevation (reach)\nW: time series of observed width (reach)\nS: time series of observed slope (reach)\nS0: bed slope\nhbf: bankfull water surface elevation\nwbf: bankfull width\nr: channel shape parameter\nz: bed elevation\n\n\n\n\n\n","category":"function"},{"location":"api/#Sad.estimate","page":"API","title":"Sad.estimate","text":"estimate(x, H, W, Qp, np, rp, zp, nens)\n\nEstimate discharge and flow parameters from SWOT observations.\n\nArguments\n\nx: channel chainage\nH: time series of water surface elevation profiles\nW:time series of width profiles\nQp: discharge prior distribution\nnp: roughness coefficient prior distribution\nrp: channel shape parameter prior distribution\nzp: downstream bed elevation prior distribution\nnens: ensemble size\nnsamples: sample size for rejection sampling\n\n\n\n\n\n","category":"function"},{"location":"algorithm/#Algorithm","page":"Algorithm","title":"Algorithm","text":"","category":"section"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"The SAD algorithm operates on the set of SWOT observables (i.e., WSE, width, and slope) and derives an estimate of river discharge and its associated uncertainty using a data assimilation scheme. The assimilation scheme involves the \"first-guess\" estimation of hydraulic variables by combining a forward model with a set of prior probability distributions before assimilating the SWOT observations. The priors are either acquired from the SWORD a-priori database or are derived using a data-driven approach (rejection sampling).","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"(Image: algorithm)","category":"page"},{"location":"algorithm/#Hydraulic-model","page":"Algorithm","title":"Hydraulic model","text":"","category":"section"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"The forward model in the assimilation scheme is based on the Gradually Varied Flow (GVF) equations, which describe the steady-state, non-uniform flow in river channels with gradual variations in water depth and velocity. The general form of the GVF equation is ","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"dfracdYdx = dfracS_0 - S_f1 - textFr^2","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"where Y is the average water depth, x is the longitudinal distance, S0 is the channel bed slope, Sf is the water surface slope, and Fr is the Froude number. The latter is given by","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"Fr = dfracQWYsqrtgY","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"where W is the average flow width, Q is river discharge, and g is the gravity acceleration, while the water surface slope can be calculated as a function of discharge, channel geometry, and flow resistance.","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"Integrating the GVF equations allows for calculation of longitudinal profiles of water surface elevation along a river, and are solved using DifferentialEquations.jl. ","category":"page"},{"location":"algorithm/#Cross-sections","page":"Algorithm","title":"Cross sections","text":"","category":"section"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"Currently, rectangular and Dingman cross sections have been implemented. The latter allows for a more realistic representation of river channel cross sections and is the default in SAD.","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"(Image: dingman)","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"The variables W^ast, W, Y_m^ast, and Y_m are the bankfull width, water surface width, bankfull maximum depth, and maximum depth respectively. These are related to the average width and flow depth via a channel shape parameter r","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"Y = left( dfracrr+1 right) Y_m","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"W = W^ast left( dfracY_mY_m^ast right)^1r","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"Hydraulic geometry relations can be incorporated into the SAD GVF model through the relationships between the AHG coefficients and exponents and the channel cross-section geometric variables: ","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"W = a Q^b","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"Y = c Q^f","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"The r parameter (0  r  infty) reflects the river channel shape, with r=1 corresponding to a triangular channel. As r increases, the channel banks become steeper, and the bottom becomes flatter leading to a rectangular channel for r rightarrow infty","category":"page"},{"location":"algorithm/#Data-assimilation","page":"Algorithm","title":"Data assimilation","text":"","category":"section"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"The assimilation algorithm employed in the implementation of SAD presented here is the Local Ensemble Transform Kalman Filter (LETKF). The LETKF is a variant of the Ensemble Kalman Filter that combines a prior probability distribution of state variables (e.g., river discharge) with direct or indirect observations (in this case, water surface elevation and width) to generate an optimal estimate (i.e., analysis). The prior distribution is represented by the model error covariance, which is calculated empirically from an ensemble of unknown model states (i.e., background ensemble). The observations and their uncertainty are represented by mapping the state variables to the observations space (e.g., river discharge to water surface elevation) and an error covariance. The analysis state (both the mean and the ensemble deviations from the mean) is essentially calculated as a function of the prior model ensemble, the model and observation error covariances, and the difference between the model-predicted observations and the actual observations.","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"The estimation of river discharge from future SWOT observations can be difficult when bed elevation and/or roughness are unknown due to equifinality. One approach that can aid in the solution of such problems is regularization, wherein additional constraints are introduced in the form of penalty terms similar to the observation difference applied to the LETKF analysis calculation. In the case of river discharge estimation, additional constraints can be derived from the at-a-station hydraulic geometry relations. In particular, it can be shown that assimilating “observations” of the form W  aQ b=0, is equivalent to a form of regularization that is adding prior knowledge (in our case adherence to the AHG equations) to help solve an ill-posed inverse problem.","category":"page"},{"location":"algorithm/#Priors","page":"Algorithm","title":"Priors","text":"","category":"section"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"Ensemble assimilation methods require the definition of a prior probability distribution from which to generate the ensemble of background variables. Given that our discharge estimation approach needs to be applicable globally, the algorithm must operate on the assumption of minimal prior knowledge regarding river discharge and the various inputs to the GVF model. The inputs to the GVF model that are not directly observed include discharge, bed elevation (as well as bed slope), the roughness coefficient, and the channel shape parameter r.","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"We adapted a rejection sampling approach to derive appropriate prior distributions for these variables. Rejection sampling is a technique used to generate samples from the target distribution T using the proposal distribution P. Instead of directly sampling from T, the method generates samples from P and accepts/rejects each of those samples according to likelihood ratio ","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"dfract(x)L p(x)","category":"page"},{"location":"algorithm/","page":"Algorithm","title":"Algorithm","text":"where L is a constant (L>1) and t(x),p(x) are the density functions of T and P, respectively. In our case, the target distribution is the prior distribution of the unobserved variable (e.g., bed elevation) and the proposal distribution is an uninformative prior. Since the density function of the target distribution is unknown, we use the GVF model as a functional to transform both densities t(x) and p(x) to correspond to density functions of water surface elevation instead of the target variable. The probability density function of WSE can be estimated from the observations, thus allowing us to calculate the likelihood ratio and accept/reject the proposed target-variable value for its prior distribution.","category":"page"},{"location":"installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"The latest version of Sad.jl can be installed with Julia's built-in package manager. You can add the package and then run instantiate to install dependencies","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"import Pkg\nPkg.add(\"Sad\")\nPkg.instantiate()","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Similarly the package can be updated ","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Pkg.update(\"Sad\")","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"warn: Use Julia 1.6 or newer\nThe Sad.jl code has only been tested with Julia v1.6 and above, so please update your installation accordingly.","category":"page"},{"location":"use_cases/#Use-cases","page":"Use cases","title":"Use cases","text":"","category":"section"},{"location":"use_cases/#Pepsi-1-experiment","page":"Use cases","title":"Pepsi-1 experiment","text":"","category":"section"},{"location":"use_cases/","page":"Use cases","title":"Use cases","text":"We start by loading the data, including river node information and synthetic SWOT observations","category":"page"},{"location":"use_cases/","page":"Use cases","title":"Use cases","text":"using NCDatasets, Statistics, Distributions, Sad\nf = Dataset(\"../../data/pepsi1/Po.nc\")\ng = NCDatasets.group(f, \"XS_Timeseries\")\nqwbm = NCDatasets.group(f, \"River_Info\")[\"QWBM\"][1]\nx = (g[\"X\"][:][end] .- g[\"X\"][:])[end:-1:1, 1]\nQ = g[\"Q\"][:][end:-1:1, :]\nH = g[\"H\"][:][end:-1:1, :]\nW = g[\"W\"][:][end:-1:1, :]","category":"page"},{"location":"use_cases/","page":"Use cases","title":"Use cases","text":"Bankfull width and water surface elevation can be guessed as the maximum from the observed time series, while an initial estimate of bed slope can be obtained from the mean of surface water slope","category":"page"},{"location":"use_cases/","page":"Use cases","title":"Use cases","text":"hbf = maximum(H, dims=2)[:, 1]\nwbf = maximum(W, dims=2)[:, 1]\nS0 = mean(diff(H, dims=1) ./ diff(x), dims=2)[:, 1]\nS0 = [S0[1]; S0]","category":"page"},{"location":"use_cases/","page":"Use cases","title":"Use cases","text":"Then we can derive the prior distributions using rejection sampling from uninformative priors","category":"page"},{"location":"use_cases/","page":"Use cases","title":"Use cases","text":"Qp0, np0, rp0, zp0 = Sad.priors(qwbm, minimum(H[1, :]), Sad.sinuous)\nQp, np, rp, zp = Sad.rejection_sampling(Qp0, np0, rp0, zp0, x, H, S0, mean(H[1, :]), wbf, hbf, 1000);","category":"page"},{"location":"use_cases/","page":"Use cases","title":"Use cases","text":"The ensemble of discharge, roughness coefficient, channel shape parameter and bed elevation can now be generated","category":"page"},{"location":"use_cases/","page":"Use cases","title":"Use cases","text":"Qe, ne, re, ze = Sad.prior_ensemble(x, Qp, np, rp, zp, 1000);","category":"page"},{"location":"use_cases/","page":"Use cases","title":"Use cases","text":"Bed elevation and slope are then estimated by assimilating the time-average water surface elevation profile","category":"page"},{"location":"use_cases/","page":"Use cases","title":"Use cases","text":"S = diff(H, dims=1) ./ diff(x);\nSe = repeat(S', outer=3)'[:, 1:1000];\nSe = [Se[1, :]'; Se]\nSad.bathymetry!(ze, Se, Qe, ne, re, x, hbf, wbf, mean(H, dims=2)[:, 1])\nzp = truncated(Normal(mean(ze[1, :]), 1e-3), -Inf, minimum(H[1, :]))\nSa = mean(Se, dims=2)[:, 1]","category":"page"},{"location":"use_cases/","page":"Use cases","title":"Use cases","text":"and finally we assimilate the observed water surface elevation to estimate discharge and flow parameters, which include roughness coefficient and minimum cross-sectional area","category":"page"},{"location":"use_cases/","page":"Use cases","title":"Use cases","text":"Qa, na = assimilate(H, W, x, wbf, hbf, Sa, Qp, np, rp, zp, nens)","category":"page"},{"location":"use_cases/","page":"Use cases","title":"Use cases","text":"(Image: po)","category":"page"},{"location":"use_cases/#Confluence","page":"Use cases","title":"Confluence","text":"","category":"section"},{"location":"#Sad.jl-Documentation","page":"Home","title":"Sad.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"installation.md\", \"algorithm.md\", \"use_cases.md\", \"api.md\"]","category":"page"},{"location":"","page":"Home","title":"Home","text":"The Surface Water Ocean Topography (SWOT) satellite mission that launched in December 2022 will offer a unique opportunity to map river discharge at an unprecedented spatial resolution globally from observations of water surface elevation, width, and slope. Because river discharge will not be directly observed from SWOT, a number of algorithms of varying complexity have been developed to estimate discharge from SWOT observables. The SAD algorithm operates on the set of SWOT observables (i.e., WSE, width, and slope) and derives an estimate of river discharge and its associated uncertainty using a data assimilation scheme. The assimilation scheme involves the “first‐guess” estimation of hydraulic variables by combining a forward model with a set of prior probability distributions before assimilating the SWOT observations. The priors are estimated with a sampling approach, as these data‐driven methods have shown promise in many fields. The objective of the algorithm is to estimate discharge at each river reach when SWOT observations become available.","category":"page"}]
}
