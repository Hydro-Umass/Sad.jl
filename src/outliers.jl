function min_slope_regression(x::Vector{Float64}, H::Vector{Union{Missing, Float64}}, min_slope::Float64)
    valid = findall(!ismissing, H)
    if isempty(valid)
        return copy(H)
    end
    z = H[valid] .- min_slope .* x[valid]
    pava_fit = Vector{Float64}(undef, length(z))
    weights = Vector{Int}(undef, length(z))
    current_block_idx = 0
    for i=1:length(z)
        current_block_idx += 1
        pava_fit[current_block_idx] = z[i]
        weights[current_block_idx] = 1
        while current_block_idx > 1 && pava_fit[current_block_idx] < pava_fit[current_block_idx - 1]
            prev_weight = weights[current_block_idx - 1]
            curr_weight = weights[current_block_idx]
            new_val = (pava_fit[current_block_idx - 1] * prev_weight + pava_fit[current_block_idx] * curr_weight) / (prev_weight + curr_weight)
            pava_fit[current_block_idx - 1] = new_val
            weights[current_block_idx - 1] += curr_weight
            current_block_idx -= 1
        end
    end
    pos = 1
    z_fit = Vector{Float64}(undef, length(z))
    for i in 1:current_block_idx
        block_size = weights[i]
        z_fit[pos:(pos + block_size - 1)] .= pava_fit[i]
        pos += block_size
    end
    h_fit = z_fit .+ min_slope .* x[valid]
    Hs = similar(H)
    Hs[valid] = h_fit
    return Hs
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
