using DifferentialEquations

"""
    froude(Q::Float64, xs::CrossSection)

Calculate the Froude number.

# Arguments
- `Q`: discharge
- `xs`: cross section

"""
function froude(Q::Float64, xs::CrossSection)
    g = 9.806
    Fr = Q / (depth(xs) * width(xs) * sqrt(g * depth(xs)))
    Fr
end

"""
    dydx(y, p, x)

Ordinary differential equation describing the Gradually-Varied-Flow model.

"""
function dydx(y, p, x)
    i = trunc(Int, ceil(x / (p[3] / length(p[1]))))
    i = i > 0 ? i : 1
    # TODO: Change this dynamically depending on whether the flow is subcritical or supercritical
    pm = -1 # integrate upstream
    Q = p[1]
    xs = p[2][i]
    xs.Ym = y * (xs.r + 1) / xs.r
    S0 = xs.S0
    Sf = xs.n^2 * Q^2 / (width(xs)^2 * depth(xs)^(10/3))
    Fr = froude(Q, xs)
    pm * (S0 - Sf) / (1 - Fr^2)
end

"""
Calculate a water surface profile by solving the Gradually-Varied-Flow equation.

- `Q`: discharge
- `ybc`: downstream boundary condition for depth
- `S0`: bed slope for reach
- `n`: roughness coefficient
- `x`: downstream distance for each cross section
- `wbf`: bankfull width for each cross section
- `ybf`: bankfull depth for each cross section
- `r`: channel geometry coefficient for each cross section

"""
function gvf(Q::Float64, ybc::Float64, S0::Vector{Float64}, n::Float64,
             x::Vector{Float64}, wbf::Vector{Float64}, ybf::Vector{Float64}, r::Vector{Float64})
    c = [Dingman(wbf[i], ybf[i], ybc, r[i], S0[i], n) for i in 1:length(x)]
    prob = DiscreteProblem(dydx, ybc, (x[1], x[end]), (Q, c, x[end]))
    h = try
        sol = solve(prob, Tsit5(), abstol=1e-2, saveat=x)
        sol.u
    catch
        zeros(length(x)) .- 9999.
            end
    h
end

"""
    gvf_ensemble!(H, S, x, hbf, wbf, Qe, ne, re, ze)

Generate an ensemble of water height profiles from Gradually-Varied-Flow simulations
and associated profiles of bed elevation.

# Arguments

- `H`: water surface elevation
- `S`: bed slope
- `x`: channel chainage
- `hbf`: bankfull water surface elevation
- `wbf`: bankfull width
- `Qe`: ensemble discharge
- `ne`: ensemble roughness coefficient
- `re`: ensemble channel shape parameter
- `ze`: ensemble bed elevation profiles

"""
function gvf_ensemble!(hbc::Float64, S, x::Vector{Float64}, hbf::Vector{Float64}, wbf::Vector{Float64}, Qe::Vector{Float64}, ne::Vector{Float64}, re::Vector{Float64}, ze)
    nens = length(Qe)
    if ndims(S) > 1
        Se = S
    else
        Se = repeat(S', outer=nens)'
    end
    for j in 2:length(x)
        ze[j, :] = ze[j-1, :] .+ Se[j, :] .* (x[j] - x[j-1]);
    end
    ybf = hbf .- ze
    he = zeros(length(x), nens)
    for i in 1:nens
        he[:, i] = gvf(Qe[i], (hbc-ze[1, i])*re[i]/(re[i] + 1), Se[:, i], ne[i], x, wbf, ybf[:, i], [re[i] for _ in 1:length(x)])
    end
    he
end

"""
    interpolate_hbc(x, H, S)

Interpolate water surface elevation boundary condition by assuming uniform flow, i.e., energy slope is equal to bed slope.

# Arguments

- `x`: channel chainage
- `H`: water surface elevation profile
- `S`: bed slope

"""
function interpolate_hbc(x::Vector{Float64}, H::Vector{FloatM}, S::Vector{Float64})
    j = minimum(findall(.!ismissing.(H)))
    h = zeros(j)
    h[1] = H[j]
    for k =2:j
        h[k] = h[k-1] - S[j+2-k] * (x[j+2-k] - x[j+1-k])
    end
    h[end]
end
