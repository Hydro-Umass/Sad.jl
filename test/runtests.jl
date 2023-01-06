using Sad
using Test
using NCDatasets, Statistics, Distributions

@testset "cross-sections" begin
    r_xs = Sad.Rectangular(100., 5., 2., 0.00001, 0.03)
    d_xs = Sad.Dingman(100., 5., 2., 2.5, 0.00001, 0.03)
    @test width(r_xs) == 100.0
    @test depth(r_xs) == 2.0
    @test width(d_xs) ≈ 69.314 atol=0.001
    @test depth(d_xs) ≈ 1.428 atol=0.001
end

# read sample data
f = Dataset("testdata.nc")
g = NCDatasets.group(f, "XS_Timeseries");
qwbm = NCDatasets.group(f, "River_Info")["QWBM"][1]
x = (g["X"][:][end] .- g["X"][:])[end:-1:1, 1]
Q = g["Q"][:][end:-1:1, :]
H = convert(Matrix{Sad.FloatM}, g["H"][:][end:-1:1, :])
W = convert(Matrix{Sad.FloatM}, g["W"][:][end:-1:1, :])
wbf = maximum.(skipmissing.(eachrow(W)))
hbf = maximum.(skipmissing.(eachrow(H)))
S = diff(H, dims=1) ./ diff(x)
S = [S[1, :]'; S]
S0 = mean(S, dims=2)[:, 1]
S = convert(Matrix{Sad.FloatM}, S)
x, H, W, S = Sad.drop_unobserved(x, H, W, S)

# perform estimation
nens = 100
nsamples = 1000
Qp0, np0, rp0, zp0 = Sad.priors(qwbm, minimum(H[1, :]), Sad.sinuous)
Qp, np, rp, zp = Sad.rejection_sampling(Qp0, np0, rp0, zp0, x, H, S0, wbf, hbf, nens, nsamples)
Qe, ne, re, ze = Sad.prior_ensemble(x, Qp, np, rp, zp, nens)

@testset "priors" begin
    @test all([p isa Distribution for p in [Qp0, np0, rp0, zp0]])
    @test all([p isa Distribution for p in [Qp, np, rp, zp]])
    @test all([mean(v) isa Real for v in [Qe, ne, re, ze]])
end

Se = Sad.bathymetry!(ze, S, S0, Qe, ne, re, x, hbf, wbf, mean.(skipmissing.(eachrow(H))))
zp = truncated(Normal(mean(ze[1, :]), 1e-3), -Inf, minimum(H[1, :]))
Sa = mean(Se, dims=2)[:, 1]

@testset "bathymetry" begin
    @test mean(ze[1, :]) isa Real
    @test mean(Sa) isa Real
end

Qa, Qu, na = Sad.assimilate(H, W, S, x, wbf, hbf, Sa, Qp, np, rp, zp, nens, nothing)
za = zeros(length(x))
za[1] = mean(zp)
for i=2:length(x) za[i] = za[i-1] + Sa[i] * (x[i] - x[i-1]) end
dA = Sad.calc_dA(Qa, na, W, S)
Wm = mean.(skipmissing.(eachcol(W)))
Sm = mean.(skipmissing.(eachcol(S)))
A0, n = Sad.flow_parameters(Qa, na, convert(Vector{Sad.FloatM}, Wm), convert(Vector{Sad.FloatM}, Sm), dA)

nse(o, m) = 1 - sum((o .- m).^2) / sum((o .- mean(o)).^2)

@testset "estimation" begin
    @test mean(Qa) > 0
    @test mean(Qu) > 0
    @test mean(Qu) < mean(Qa)
    @test mean(na) > 0.01 && mean(na) < 0.07
    @test A0 > 0
    @test n > 0.01 && n < 0.07
    @test nse(Q[1, :], Qa) > 0.9
end
