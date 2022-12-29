using Sad: CrossSection


"""
Rectangular river channel cross section.

- `Wb`: bankfull width
- `Yb`: bankfull depth
- `Ym`: flow depth
- `S0`: bed slope
- `n`: Manning roughness coefficient

See also: [`Dingman`](@ref)
"""
mutable struct Rectangular <: CrossSection
    Wb:: Real
    Yb:: Real
    Ym:: Real
    S0:: Real
    n:: Real
end

"""
Dingman river channel cross section.

Uses generalized expressions for cross-section geometry based on at-a-station hydraulic geometry relations.

- `Wb`: bankfull width
- `Yb`: bankfull depth
- `Ym`: flow depth
- `r`: geometry exponent
- `S0`: bed slope
- `n`: Manning roughness coefficient

See also: [`Rectangular`](@ref)
"""
mutable struct Dingman <: CrossSection
    Wb:: Real
    Yb:: Real
    Ym:: Real
    r:: Real
    S0:: Real
    n:: Real
end

"""
    width(xs)

Calculate water surface width of cross section.

# Arguments
- `xs`: cross section

"""
function width(xs::Dingman)
    xs.Wb * (xs.Ym / xs.Yb)^(1 / xs.r)
end

function width(xs::Rectangular)
    xs.Wb
end

"""
    depth(xs)

Calculate average depth of cross section.

# Arguments
- `xs`: cross section

"""
function depth(xs::Dingman)
    (xs.r / (xs.r + 1)) * xs.Ym
end

function depth(xs::Rectangular)
    xs.Ym
end

export width, depth
