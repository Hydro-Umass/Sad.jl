using LinearAlgebra
using Statistics

"""
    letkf(A::Matrix{Float64}, d::Vector{Float64}, HA::Matrix{Float64}, E::Matrix{Float64}, ρ::Float64, diagR::Bool)

Apply the Local Ensemble Transform Kalman Filter algorithm.

# Arguments

- `A`: state variable ensemble matrix
- `d`: observation vector
- `HA`: model-predicted observation ensemble matrix
- `E`: observation error ensemble matrix
- `ρ`: covariance inflation parameter (optional, default value of 1.05)
- `diagR`: force observation error covariance to be diagonal, i.e. no spatial correlation (optional, default is false)

"""
function letkf(A::Matrix{Float64}, d::Vector{Float64}, HA::Matrix{Float64}, E::Matrix{Float64}; ρ::Float64=1.05, diagR::Bool=false)
    Aa = zeros(size(A))
    _, nens = size(A)
    nobs = length(d)
    if size(E) == (nobs, nobs)
        R = E
    else
        if diagR
            R = diagm(0 => diag((1 / (nens - 1)) * E * E'))
        else
            R = (1 / (nens - 1)) * E * E'
        end
    end
    Y = HA .- mean(HA, dims=2)
    X = A .- mean(A, dims=2)
    # FIXME restore localization for the analysis step
    C = Y' * pinv(R)
    P = pinv((nens - 1) / ρ * I + C * Y)
    W = real(sqrt((nens - 1) * P))
    w = P * C * (d .- mean(HA, dims=2))
    W = W .+ w
    Xᵃ = X * W .+ mean(A, dims=2)
    Aa = Xᵃ
    Aa
end
