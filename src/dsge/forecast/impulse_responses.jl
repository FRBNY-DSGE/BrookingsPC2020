import DSGE: impulse_responses

function impulse_responses(m::AbstractDSGEModel, paras::Matrix{S},
                           input_type::Symbol, method::Symbol,
                           output_vars::Vector{Symbol} =
                           [:irfstates, :irfobs, :irfpseudo]; parallel::Bool = false,
                           permute_mat::Matrix{S} = Matrix{Float64}(undef,0,0),
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           n_obs_var::Int = 0,
                           flip_shocks::Bool = false,
                           cholesky_obs_shock::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           compute_meansbands::Bool = false,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           do_cond_obs_shocks::Bool = false,
                           do_rev_transform::Bool = false,
                           correct_integ_series::Bool = true,
                           verbose::Symbol = :high) where {S<:Real}
    if n_obs_var <= 0
        error("To use method $method, user must specify the index of" *
              " the target observable with keyword n_obs_var.")
    end

    # Set up computation method
    mapfcn = parallel ? pmap : map
    h = impulse_response_horizons(m)
    fcn = if method == :max_business_cycle_variance || method == :maxBC
        function _maxBCirf(model, para)
            DSGE.update!(model, para)
            system = compute_system(model)
            states, obs, pseudo = impulse_responses(system, h, frequency_band,
                                                    n_obs_var, flip_shocks =
                                                    flip_shocks)
            if do_rev_transform
                for (k,v) in m1.observable_mappings
                    irf_trans = DSGE.get_irf_transform(v.rev_transform)
                    obs[m1.observables[k],:] = irf_trans(obs[m1.observables[k],:])
                end
            end
            return states, obs, pseudo
        end
    elseif method == :cholesky
        if isempty(permute_mat)
            @warn "Permutation matrix permute_mat is empty. Defaulting to identity."
            permute_mat = Matrix{S}(I, n_observables(m), n_observables(m))
        end
        function _choleskyirf(model, para)
            DSGE.update!(model, para)
            system = compute_system(model)

            obs_shock = zeros(n_observables(model))
            obs_shock[n_obs_var] = 1.
            _, obs, _ = impulse_responses(model, system, h, permute_mat,
                                                    obs_shock, flip_shocks =
                                                    flip_shocks,
                                                    cholesky_obs_shock = cholesky_obs_shock)
            obs_shock[n_obs_var] /= obs[n_obs_var,1]
            states, obs, pseudo = impulse_responses(model, system, h, permute_mat,
                                                    obs_shock, flip_shocks =
                                                    flip_shocks,
                                                    cholesky_obs_shock = cholesky_obs_shock)

            if do_rev_transform
                for (k,v) in m1.observable_mappings
                    irf_trans = DSGE.get_irf_transform(v.rev_transform)
                    obs[m1.observables[k],:] = irf_trans(obs[m1.observables[k],:])
                end
            end
            return states, obs, pseudo
        end
    elseif method == :regime_switching_structural_irfs
        function _struct_irfs(model, para)
            DSGE.update!(model, para)
            regime_system = compute_system(model; regime_switching = true)
            n_regs = n_regimes(regime_system)
            vec_irfstates = Vector{Array{S}}(undef,n_regs)
            vec_irfobs    = Vector{Array{S}}(undef,n_regs)
            vec_irfpseudo = Vector{Array{S}}(undef,n_regs)
            for i = 1:n_regimes(regime_system)
                vec_irfstates[i], vec_irfobs[i], vec_irfpseudo[i] =
                    impulse_responses(m, System(regime_system, i))
            end
            return vec_irfstates, vec_irfobs, vec_irfpseudo
        end
    end

    # Compute IRFs
    paras = mapslices(x -> [vec(x)], paras, dims = 2)
    irf_output = mapfcn(para -> fcn(m, para), paras)

    # Reformat output
    states = map(x -> x[1], irf_output)
    obs    = map(x -> x[2], irf_output)
    pseudo = map(x -> x[3], irf_output)

    if compute_meansbands
        # Set up metadata and output from IRFs computation
        mb_output_vars = Vector{Matrix{Matrix{S}}}(undef,0)
        class_vars = sort(map(x -> get_class(x), output_vars))
        if :obs in class_vars
            push!(mb_output_vars, obs)
        end
        if :pseudo in class_vars
            push!(mb_output_vars, pseudo)
        end
        if :states in class_vars
            push!(mb_output_vars, states)
        end

        if method != :regime_switching_structural_irfs
            mb_vec = Vector{MeansBands}(undef,0)
            for (mb_i, output_var, class) in zip(1:length(class_vars), mb_output_vars, class_vars)
                metadata = Dict{Symbol,Any}()
                metadata[:para] = input_type
                metadata[:cond_type] = :none
                metadata[:product] = :irf
                metadata[:class] = class
                metadata[:date_inds] = OrderedDict()

                # Set up for loop over variable names
                means = DataFrame()
                bands = Dict{Symbol,DataFrame}()
                names = if class == :states
                    keys(m.endogenous_states)
                elseif class == :obs
                    keys(m.observables)
                elseif class == :pseudo
                    keys(m.pseudo_observables)
                end
                metadata[:indices] = OrderedDict{Symbol,Int}(name => name_i
                                                             for (name_i, name) in enumerate(names))

                # Means and Bands for each variable in a class
                for (name_i,name) in enumerate(names)
                    # Note: obs is Vector{nobs x nperiods} -> for each observable,
                    # so we want to select the IRF of a specific obs, i.e. map(x -> x[obs_index,:]).
                    # This creates a nperiod x ndraws matrix, which we want to transpose
                    # to get a ndraws x nperiod matrix
                    single_var = Matrix(reduce(hcat, map(x -> x[name_i,:], output_var))')
                    if name in [:laborproductivity, :Wages] && correct_integ_series
                        zt_tmp = Matrix(reduce(hcat, map(x -> x[findfirst(:z_t .== names),:],
                                                         output_var))')
                        zt_tmp = cumsum(zt_tmp, dims = 2)
                        single_var += zt_tmp
                    end
                    means[!,name] = vec(mean(single_var, dims = 1))
                    bands[name]   = find_density_bands(single_var, density_bands;
                                                       minimize = minimize)
                end
                push!(mb_vec, MeansBands(metadata, means, bands))

                # Save MeansBands
                tail = method == :cholesky ? :_cholesky : :_maxBC
                fp = get_meansbands_output_file(m, input_type, :none,
                                                Symbol(:irf, class, tail),
                                                forecast_string = forecast_string,
                                                do_cond_obs_shocks = do_cond_obs_shocks)
                dirpath = dirname(fp)
                isdir(dirpath) || mkpath(dirpath)
                JLD2.jldopen(fp, true, true, true, IOStream) do file
                    write(file, "mb", mb_vec[mb_i])
                end
                println(verbose, :high, "  " * "wrote " * basename(fp))
            end
        else
        end
        return mb_vec
    else
        # Reshape to be nobs x nperiod x ndraw
        if method == :regime_switching_structural_irfs
            n_regs = length(states[1]) # since states is a Vector of Vector{Array}s,
                                       # where the Vector{Array}s is a vector of
                                       # impulse response matrices.
            states = cat(map(x -> [x[i]' for i = 1:n_regs], states)..., dims = 3)
            obs    = cat(map(x -> [x[i]' for i = 1:n_regs], obs)..., dims = 3)
            pseudo = cat(map(x -> [x[i]' for i = 1:n_regs], pseudo)..., dims = 3)
        else
            states = cat(map(x -> x', states)..., dims = 3)
            obs    = cat(map(x -> x', obs)..., dims = 3)
            pseudo = cat(map(x -> x', pseudo)..., dims = 3)
        end
        return states, obs, pseudo
    end
end

function long_run_cholesky_impulse_responses(m::AbstractRepModel,
                                             system::System{S}, horizon::Int64;
                                             flip_shocks::Bool = false,
                                             get_shocks::Bool = false,
                                             use_pinv::Bool = false) where {S<:Real}
    # Set up IRFs
    nshocks = size(system[:RRR], 2)
    nstates = size(system[:TTT], 1)
    nobs    = size(system[:ZZ], 1)

    system = DSGE.zero_system_constants(system) # Set constant system matrices to 0
    s₀ = zeros(S, nstates)

    # Back out structural shocks corresponding to long-run Cholesky shock to observables
    lr_obs_std = system[:ZZ] * inv(Matrix{S}(I, nstates, nstates) - system[:TTT]) * system[:RRR]
    lr_obs_cov = lr_obs_std * system[:QQ] * lr_obs_std'
    cholmat = if nobs != nshocks || use_pinv
        pinv(lr_obs_std * sqrt.(system[:QQ])) *
            cholesky((lr_obs_cov + lr_obs_cov') ./ 2).U'
    else
        inv(lr_obs_std * sqrt.(system[:QQ])) *
            cholesky((lr_obs_cov + lr_obs_cov') ./ 2).U'
    end
    struct_shock = cholmat[:,1]

    # Run IRF
    shocks_augmented = zeros(S, nshocks, horizon)
    shocks_augmented[:,1] = flip_shocks ? struct_shock : -struct_shock # negative shock by default
    states, obs, pseudo = forecast(system, s₀, shocks_augmented)

    if get_shocks
        return states, obs, pseudo, struct_shock, lr_obs_cov, lr_obs_std
    else
        return states, obs, pseudo
    end
end
