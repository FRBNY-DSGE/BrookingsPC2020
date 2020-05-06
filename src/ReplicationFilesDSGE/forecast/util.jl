# This function zeros out shocks not specified by shocks
# and updates the AbstractDSGEModel object accordingly
function zero_other_shocks!(m::AbstractDSGEModel{S}, shocks::Vector{Symbol}) where {S<:Real}
    si = map(k -> m.exogenous_shocks[k], shocks) # shock indices
    zero_shocks = setdiff(collect(keys(m.exogenous_shocks)), shocks)
    for shock in zero_shocks
        σ_name = Symbol(:σ_, replace(string(shock), "_sh" => ""))
        if σ_name == :σ_zp
            σ_name = :σ_z_p
        elseif σ_name == :σ_rm
            σ_name = :σ_r_m
        elseif string(shock)[end-1:end] != "sh"
            σ_name = Symbol(:σ_r_m, string(shock)[end])
        end
        m[σ_name].value = zero(S)
    end
end

# This function computes the corresponding transition and measurement
# equations specified by the desired observables and shocks
# using the original ZZ and ZZ_pseudo
function update_system!(m::AbstractDSGEModel{S}, system::System,
                        observables::Vector{Symbol}, shocks::Vector{Symbol};
                        zero_DD::Bool = false) where {S<:Real}
    # Set up indices
    oid = m.observables # observables indices dictionary
    pid = m.pseudo_observables # pseudo observables indices dictionary

    # Compute new ZZ and DD matrices
    Zout = zeros(S, length(observables), n_states_augmented(m))
    Dout = zeros(S, length(observables))
    Eout = zeros(S, length(observables), length(observables))
    for (i,obs) in enumerate(observables)
        Zout[i,:], Dout[i] = if haskey(oid, obs)
            system[:ZZ][oid[obs], :], zero_DD ?
                zero(S) : system[:DD][oid[obs]]
        elseif haskey(pid, obs)
            system[:ZZ_pseudo][pid[obs], :], zero_DD ? zero(S) : system[:DD_pseudo][pid[obs]]
        else
            error("Observable/PseudoObservable $obs cannot be found in the DSGE model $m")
        end
    end

    # Find shocks to keep
    shock_inds = map(k -> m.exogenous_shocks[k], shocks)

    # Update system
    system.measurement.ZZ = Zout
    system.measurement.DD = Dout
    system.measurement.QQ = system[:QQ][shock_inds, shock_inds]
    system.measurement.EE = Eout
    system.transition.RRR = system[:RRR][:, shock_inds]
end
