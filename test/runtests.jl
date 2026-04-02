using Sad
using Test
using NCDatasets, Statistics, Distributions

# read sample data
f = Dataset("testdata.nc")
g = f.group["XS_Timeseries"]
qwbm = f.group["River_Info"]["QWBM"][1]
x = (g["X"][:][end] .- g["X"][:])[end:-1:1, 1]
Q = g["Q"][end:-1:1, :]
H = convert(Matrix{Sad.FloatM}, g["H"][end:-1:1, :])
W = convert(Matrix{Sad.FloatM}, g["W"][end:-1:1, :])
S = diff(H, dims=1) ./ diff(x)
S = convert(Matrix{Sad.FloatM}, [S[1, :]'; S])

reach = Sad.preprocess(x, H, W, S)
p = Sad.priors(qwbm, reach.hmin, Sad.sinuous)
res = Sad.infer(p, reach)

@testset "inference" begin
    @test
W[W .== -9999.0] .= missing
S = diff(H, dims=1) ./ diff(x)
S = [S[1, :]'; S]
x, H, W, S = Sad.drop_unobserved(x, H, W, S)
Qp0, np0, rp0, zp0 = Sad.priors(qwbm, minimum(skipmissing(H[1, :])), Sad.braided)
Qa, Qu, A0, n = Sad.estimate(x, H, W, S, nothing, Qp0, np0, rp0, zp0, nens, nsamples, nothing, nothing, nothing)

@testset "missing data" begin
    @test length(findall(ismissing.(Qa))) == 315
    @test A0 > 0
    @test n > 0.01 && n < 0.07
end
