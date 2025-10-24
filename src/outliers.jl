using DataInterpolations

function monotonic_interpolate_min_slope(x::Vector{Float64}, y::Vector{FloatM}, xq::Vector{Float64}; min_slope::Float64=1e-4)
    idx = sortperm(x)
    x, y = x[idx], y[idx]
    good = .!ismissing.(y)
    if sum(good) > 2
        itp = PCHIPInterpolation(y[good], x[good]; extrapolation=ExtrapolationType.Linear)
        yq = convert(Vector{FloatM}, itp.(xq))
    else
        yq = copy(y)
    end
    j = findall(.!ismissing.(y))
    for i in 2:length(j)
        dx = xq[j[i]] - xq[j[i-1]]
        required = yq[j[i-1]] + min_slope * dx
        if yq[j[i]] <= required
            yq[j[i]] = required
        end
    end
    yq[ismissing.(y)] .= missing
    return yq
end

function calc_slope(x::Vector{Float64}, H::Vector{FloatM}, min_slope::Float64)
    S = diff(H) ./ diff(x)
    S = vcat(S, S[end])
    # ensure that the slope vector has the correct type in case all slopes are missing
    S = convert(Vector{FloatM}, S)
    invalid = findall(ismissing.(S) .& .!ismissing.(H))
    if length(invalid) > 0
        j = findall(.!ismissing.(H))
        if length(j) > 1
            S[j[1:end-1]] = diff(H[j]) ./ diff(x[j])
            S[j[end]] = S[j[end-1]]
        else
            S[j] .= min_slope
        end
    end
    return S
end
