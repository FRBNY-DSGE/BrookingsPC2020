function ait()
    AltPolicy(:zero_pi_star, zero_pi_star_eqcond, zero_pi_star_solve,
              forecast_init = zero_pi_star_forecast_init,
              color = RGB(0., 0., 0.5430)) # dark blue
end

"""
```
zero_pi_star_eqcond(m::AbstractDSGEModel)
```

Solves for the transition equation of `m` under a price level
targeting rule (implemented by adding a price-gap state)
"""
function eqcond_zero_pi_star(m::AbstractDSGEModel)

    # get the old indices
    old_states = sort!(collect(values(m.endogenous_states)))
    old_eqs    = sort!(collect(values(m.equilibrium_conditions)))

    # Get equilibrium condition matrices
    Γ0_noaltpol, Γ1_noaltpol, C_noaltpol, Ψ_noaltpol, Π_noaltpol  = eqcond(m)

    # check to make sure this state hasn't already been added, then:
    # 1. add the pgap state to the model
    # 2. increment the indices of every augmented state

    nstates = n_states(m)

    # fill in new Γ0, Γ1, C, Ψ, Π
    Γ0 = zeros(Float64, nstates, nstates)
    Γ0[old_eqs, old_states] = Γ0_noaltpol

    Γ1 = zeros(Float64, nstates, nstates)
    Γ1[old_eqs, old_states] = Γ1_noaltpol

    C  = zeros(Float64, nstates)
    C[old_eqs, :] = C_noaltpol

    Ψ  = zeros(Float64, nstates, n_shocks_exogenous(m))
    Ψ[old_eqs, :] = Ψ_noaltpol

    Π  = zeros(Float64, nstates, n_shocks_expectational(m))
    Π[old_eqs, :] = Π_noaltpol

    # add law of motion for pgap
    eq             = m.equilibrium_conditions
    endo           = m.endogenous_states

    # This assumes that the inflation target is the model's steady state
    Γ0[eq[:eq_mp], endo[:π_star_t]]  =  0.0
    Γ0[eq[:eq_π_star], endo[:π_star_t]] = 1.0
    Γ1[eq[:eq_π_star], endo[:π_star_t]] = 0.0

    return Γ0, Γ1, C, Ψ, Π
end

"""
```
zero_pi_star_solve(m::AbstractDSGEModel)
```

Solves for the transition equation of `m` under a price level
targeting rule (implemented by adding a price-gap state)
"""
function solve_zero_pi_star(m::AbstractDSGEModel)

    # Get equilibrium condition matrices
    Γ0, Γ1, C, Ψ, Π  = eqcond_zero_pi_star(m)

    # Solve model
    TTT_gensys, CCC_gensys, RRR_gensys, eu = gensys(Γ0, Γ1, C, Ψ, Π, 1+1e-6)

    if !((eu[1] == 1) & (eu[2] == 1))
        throw(GensysError("Gensys does not give existence"))
    end
    TTT_gensys = real(TTT_gensys)
    RRR_gensys = real(RRR_gensys)
    CCC_gensys = reshape(CCC_gensys, size(CCC_gensys, 1))

    # Augment states
    TTT, RRR, CCC = DSGE.augment_states(m, TTT_gensys, RRR_gensys, CCC_gensys)

    return TTT, RRR, CCC
end

"""
```
init_zero_pi_star_forecast(m::AbstractDSGEModel, shocks::Matrix{T}, final_state::Vector{T})
```

Adjust shocks matrix and final state vector for forecasting under the ZERO_PI_STAR rule
"""
function forecast_init_zero_pi_star(m::AbstractDSGEModel, shocks::Matrix{T}, final_state::Vector{T};
                           cond_type::Symbol = :none) where {T<:AbstractFloat}
      return shocks, final_state
end
