import DSGE: impulse_responses

# This file holds the VAR versions of impulse responses.
function impulse_responses(beta::Matrix{S}, sigma::Matrix{S}, n_obs_shock::Int,
                           horizon::Int, shock_size::S = one(S);
                           method::Symbol = :cholesky,
                           flip_shocks::Bool = false,
                           beta_has_constant::Bool = true,
                           frequency_band::Tuple{S,S} =
                           (2*π/32, 2*π/6)) where {S<:Real}

    # Compute dimensions
    n = size(beta,2)
    lags = convert(Int, beta_has_constant ? (size(beta,1) - 1) / n : size(beta,1) / n)

    # Compute impact based on IRF type
    Y = zeros(lags + horizon, n)
    Y[lags + 1, :] = if method == :cholesky
        cholesky_shock(sigma, n, n_obs_shock, shock_size;
                       flip_shocks = flip_shocks)
    elseif method == :maximum_business_cycle_variance || method == :maxBC
        maxBC_shock(beta, sigma, n, n_obs_shock, shock_size, lags, frequency_band;
                    flip_shocks = flip_shocks)
    elseif method == :choleskyLR || method == :cholesky_long_run
        cholesky_long_run_shock(beta, sigma, n_obs_shock, n, lags, shock_size;
                                flip_shocks = flip_shocks)
    else
        error("IRF method $(string(method)) has not been implemented.")
    end

    # For efficiency
    if beta_has_constant
        beta = @views beta[2:end, :]
    end

    # Compute impulse response
    for t = 2:horizon
        xT = reshape(Y[lags + t - 1:-1:lags + t - lags, :]', lags * n, 1)'
        Y[lags + t, :] = xT * beta
    end

    return Y[lags + 1:end, :]
end


function cholesky_shock(sigma::Matrix{S}, n::Int, n_obs_shock::Int,
                        shock_size::S; flip_shocks::Bool = false) where {S<:Real}
    cholmat = cholesky((sigma + sigma') ./ 2).L
    vec_shock = zeros(n)
    vec_shock[n_obs_shock] = flip_shocks ? shock_size : -shock_size # negative by DSGE convention
    return (cholmat * vec_shock)'
end

function maxBC_shock(beta::Matrix{S}, sigma::Matrix{S}, n::Int, n_obs_shock::Int, shock_size::S,
                     lags::Int, frequency_band::Tuple{S,S};
                     flip_shocks::Bool = false) where {S<:Real}
    if lags * n < size(beta,1)
        beta = @views beta[2:end, :]
    end

    cholmat = cholesky((sigma + sigma') ./ 2).L
    increment = abs(frequency_band[1] - frequency_band[2]) / 200.
    V = zeros(S, n, n) # variance
    eminusif = zeros(Complex{S}, 1, 1, lags)
    for f = frequency_band[1]:increment:round(frequency_band[2], digits=10) # rounding difference from Matlab sometimes leads to one fewer loop than Matlab computes
        eminusif[1, 1, :] = exp.(-im .* f .* collect(1:lags))
        sumB = dropdims(sum(reshape(beta', n, n, lags) .*
                           repeat(eminusif, n, n, 1); dims = 3), dims = 3)
        invA = (Matrix{Complex{S}}(I, n, n) - sumB) \ cholmat
        V += reshape(real.(kron(conj(invA[n_obs_shock, :]), invA[n_obs_shock, :])), n, n) .*
            increment ./ abs(frequency_band[1] - frequency_band[2])
    end
    eigout = eigen(V)
    q = eigout.vectors[:, argmax(eigout.values)]
    q .*= sign(q[n_obs_shock])
    q .*= flip_shocks ? shock_size : -shock_size # negative by DSGE convention

    return (cholmat * q)'
end

function cholesky_long_run_shock(β::Matrix{S}, Σ::Matrix{S}, n_obs_shock::Int, n::Int,
                                 lags::Int, shock_size::S;
                                 flip_shocks::Bool = false) where {S<:Real}
    if n * lags < size(β, 1)
        β = β[2:end,:] # don't need the constant
    end

    # Compute decomposition
    B̃ = Matrix{S}(I, n, n) - dropdims(sum(reshape(β', n, n, lags), dims = 3), dims = 3)
    S̃ = B̃ \ (Σ * inv(B̃)')             # LR covariance = S̃ = B̃⁻¹ * Σ * B̃⁻¹' =>
    Γ = B̃ * cholesky((S̃ + S̃') ./ 2).L # S = B̃ \ (Σ * B̃⁻¹')

    # Compute shock
    vec_shock = zeros(n)
    vec_shock[n_obs_shock] = flip_shocks ? shock_size : -shock_size # negative by DSGE convention
    return (Γ * vec_shock)'
end
