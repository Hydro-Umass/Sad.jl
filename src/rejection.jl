using Statistics

"""
    coverage(reach, t) -> Int

Count non-missing H observations at timestep t across all SWOT nodes.
More missing nodes = less spatial information for constraining reach parameters.
"""
function coverage(reach::SWOTReach, t::Int)
    sum(.!ismissing.(reach.obs.H[:, t]))
end

"""
    select_representative_timesteps(reach; n_bins=5) -> Vector{Int}

Select timesteps that cover the range of downstream WSE by quantile binning
and picking the timestep with the most complete spatial coverage per bin.

Quantile binning ensures equal representation of low, mid, and high flow
conditions — equal-width bins underrepresent high flow where WSE varies
slowly with Q.
"""
function select_representative_timesteps(reach::SWOTReach; n_bins::Int=5)
    valid_ts = findall(reach.valid)
    isempty(valid_ts) && return Int[]

    h_ds = [reach.H[1, t] for t in valid_ts]

    if minimum(h_ds) ≈ maximum(h_ds)
        best = valid_ts[argmax(coverage.(Ref(reach), valid_ts))]
        @warn "All downstream WSE values identical — returning single timestep"
        return [best]
    end

    n_bins  = min(n_bins, length(valid_ts))
    quantiles = range(0.0, 1.0, length=n_bins+1)
    edges     = quantile(h_ds, quantiles)
    selected  = Int[]

    for i in 1:n_bins
        lo, hi = edges[i], edges[i+1]
        in_bin = if i == n_bins
            findall(h -> lo <= h <= hi, h_ds)
        else
            findall(h -> lo <= h < hi, h_ds)
        end
        isempty(in_bin) && continue
        best = in_bin[argmax(coverage.(Ref(reach), valid_ts[in_bin]))]
        push!(selected, valid_ts[best])
    end

    return selected
end

"""
    rejection_sample(priors, reach; N, ε_rel, max_attempts, n_bins)

Draw samples from `priors` and reject those whose GVF-predicted WSE
deviates from observations by more than `ε_rel` (mean relative RMSE
across representative timesteps).

Returns a (4 × N_accepted) matrix with rows [n, r, z0, α].
Q is not included — it is handled per-timestep by q_ensemble.
"""
function rejection_sample(priors::SWOTPriors, reach::SWOTReach,
                          months::Vector{Int} = fill(1, reach.nt);
                          N::Int            = 500,
                          ε_rel::Float64    = 0.5,
                          max_attempts::Int = 20_000,
                          n_bins::Int       = 5)
    rep_ts = select_representative_timesteps(reach; n_bins=n_bins)
    if isempty(rep_ts)
        @warn "No valid timesteps found for rejection sampling"
        return Matrix{Float64}(undef, 4, 0)
    end

    h_ds  = reach.H[1, reach.valid]
    H_ref = max(maximum(h_ds) - minimum(h_ds), 1.0)

    obs_per_t = map(rep_ts) do t
        good = findall(.!ismissing.(reach.obs.H[:, t]))
        (x    = collect(reach.obs.x[good]),
         H    = Float64.(reach.obs.H[good, t]),
         H_bc = reach.H[1, t])
    end

    accepted = Matrix{Float64}(undef, 4, N)
    n_acc    = 0
    n_try    = 0

    while n_acc < N && n_try < max_attempts
        n_try += 1

        n_draw = rand(priors.np)
        r_draw = rand(priors.rp)
        z0     = rand(priors.zp)
        α      = rand(priors.ap)

        rmse_all = Float64[]
        failed   = false

        for (obs, t) in zip(obs_per_t, rep_ts)
            Qp_t = monthly_q_prior(priors, months[t])
            Q_t  = rand(Qp_t)
            pred = gvf_solve(Q_t, n_draw, r_draw, α, z0, obs.H_bc, reach;
                             saveat=obs.x)
            if isnothing(pred)
                failed = true
                break
            end
            rmse = sqrt(mean((pred .- obs.H).^2))
            push!(rmse_all, rmse)
        end

        failed && continue
        mean(rmse_all) / H_ref > ε_rel && continue

        n_acc += 1
        accepted[:, n_acc] = [n_draw, r_draw, z0, α]
    end

    rate = round(100 * n_acc / n_try, digits=1)
    @info "Rejection sampling: accepted $n_acc / $n_try ($rate%) across $(length(rep_ts)) representative timesteps"
    if rate > 20.0
        @warn "Acceptance rate > 20% — consider tightening ε_rel"
    end
    if rate < 1.0 && n_acc < N
        @warn "Acceptance rate < 1% — consider loosening ε_rel"
    end

    return accepted[:, 1:n_acc]
end

"""
    q_ensemble(priors, reach_ensemble, reach, t; eps_abs, n_tries)

Find Q values consistent with the observed upstream WSE at timestep t,
conditioned on the accepted reach parameter ensemble.

Uses the most upstream valid node for conditioning — the downstream WSE
is the ODE boundary condition and carries no Q information. Samples Q
uniformly across [ppf(0.01), ppf(0.999)] to ensure the ensemble spans
the full plausible range including extreme flood events.

Returns (Q_accepted, member_idx) where member_idx indexes into
reach_ensemble columns, preserving the (Q, n, r, z0, α) pairing.
Returns (Float64[], Int[]) if no Q values are accepted.
"""
function q_ensemble(priors::SWOTPriors,
                    reach_ensemble::Matrix{Float64},
                    reach::SWOTReach,
                    t::Int,
                    month::Int = 1;
                    eps_abs::Float64 = 1.0,
                    n_tries::Int     = 10)
    empty_result = (Float64[], Int[])

    !reach.valid[t] && return empty_result

    H_col       = reach.obs.H[:, t]
    valid_nodes = findall(.!ismissing.(H_col))
    length(valid_nodes) < 2 && return empty_result

    upstream_node = valid_nodes[end]   # most upstream valid node
    H_obs_up      = Float64(H_col[upstream_node])
    x_up          = [reach.obs.x[upstream_node]]
    H_bc          = reach.H[1, t]
    N             = size(reach_ensemble, 2)

    # scale acceptance window to 10% of WSE range at the upstream node
    h_up_valid = [Float64(reach.obs.H[upstream_node, tt])
                  for tt in findall(reach.valid)
                  if !ismissing(reach.obs.H[upstream_node, tt])]
    eps_use = if length(h_up_valid) >= 2
        0.1 * max(maximum(h_up_valid) - minimum(h_up_valid), 1.0)
    else
        eps_abs
    end

    # sample Q uniformly across the monthly prior support
    Qp_t = monthly_q_prior(priors, month)
    Q_lo = quantile(Qp_t, 0.01)
    Q_hi = quantile(Qp_t, 0.999)

    accepted_Q   = Dict{Int, Float64}()
    tile_idx     = repeat(1:N, n_tries)

    for (i, member) in enumerate(tile_idx)
        haskey(accepted_Q, member) && continue

        Q_try = Q_lo + rand() * (Q_hi - Q_lo)
        n_k, r_k, z0_k, α_k = reach_ensemble[:, member]

        pred = gvf_solve(Q_try, n_k, r_k, α_k, z0_k, H_bc, reach;
                         saveat=x_up)
        isnothing(pred) && continue

        if abs(pred[1] - H_obs_up) < eps_use
            accepted_Q[member] = Q_try
        end
    end

    if isempty(accepted_Q)
        @warn "q_ensemble: no Q accepted at timestep $t"
        return empty_result
    end

    members = collect(keys(accepted_Q))
    Q_out   = [accepted_Q[m] for m in members]
    return Q_out, members
end
