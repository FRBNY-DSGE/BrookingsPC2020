function cholesky_responses(system::System{S}, horizon::Int64,
                           shock_vector::Vector{S};
                           flip_shocks::Bool = false,
                           get_shocks::Bool = false) where {S<:Real}
    # Set up IRFs
    nshocks = size(system[:RRR], 2)
    nstates = size(system[:TTT], 1)
    nobs    = size(system[:ZZ], 1)
    npseudo = size(system[:ZZ_pseudo], 1)

    system = DSGE.zero_system_constants(system) # Set constant system matrices to 0
    s₀ = zeros(S, nstates)

    # Compute Cholesky shock
    obs_cov_mat = (system[:ZZ] * system[:RRR]) * system[:QQ] *
        (system[:ZZ] * system[:RRR])'
    cholmat1 = cholesky((obs_cov_mat + obs_cov_mat') ./ 2).U'
    cholmat2 = pinv(system[:ZZ] * system[:RRR] * sqrt.(system[:QQ])) * cholmat1

    deviation1 = cholmat1 * shock_vector
    deviation2 = system[:RRR] * sqrt.(system[:QQ]) * cholmat2 * shock_vector

    return deviation1, deviation2
end
