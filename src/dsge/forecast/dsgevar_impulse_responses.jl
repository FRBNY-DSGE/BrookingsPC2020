import DSGE: impulse_responses

function impulse_responses(m::AbstractDSGEModel, paras::Matrix{S},
                           input_type::Symbol, method::Symbol,
                           lags::Int, observables::Vector{Symbol},
                           exogenous_shocks::Vector{Symbol},
                           n_obs_var::Int;
                           parallel::Bool = false,
                           regime_switching::Bool = false,
                           n_regimes::Int = 2,
                           regime::Int = 1,
                           frequency_band::Tuple{S,S} = (2*π/32, 2*π/6),
                           flip_shocks::Bool = false,
                           density_bands::Vector{Float64} = [.5, .6, .7, .8, .9],
                           compute_meansbands::Bool = false,
                           annualize_inflation::Bool = true,
                           minimize::Bool = true,
                           forecast_string::String = "",
                           addl_text::String = "",
                           correct_integ_series::Bool = false,
                           compute_inflation_expectation_error::Bool = true,
                           verbose::Symbol = :high) where {S<:Real}
    if n_obs_var <= 0
        error("To use method $method, user must specify the index of" *
              " the target observable with keyword n_obs_var.")
    end

    if correct_integ_series &  !(:z_t in observables)
        error("Productivity process z_t must be an observable to correct integrated time series.")
    end

    # Set up computation method
    mapfcn = parallel ? pmap : map
    h = impulse_response_horizons(m)

    # Compute VAR coefficients implied by DSGE
    paras = mapslices(x -> [vec(x)], paras, dims = 2)
    var_output =
        mapfcn(para -> DSGE_to_VAR!(m, para, observables, shocks, lags;
                                    regime_switching = regime_switching,
                                    n_regimes = n_regimes,
                                    regime = regime), paras)

    # Reformat output
    β_draws = map(x -> x[1], var_output)
    Σ_draws = map(x -> x[2], var_output)

    # Compute IRFs
    irf_output =
        mapfcn((β, Σ) ->
               impulse_responses(β, Σ, n_obs_var, h;
                                 method = method,
                                 beta_has_constant = false,
                                 flip_shocks = flip_shocks),
               β_draws, Σ_draws)

    if compute_meansbands
        # Set up metadata and output from IRFs computation
        metadata = Dict{Symbol,Any}()
        metadata[:para] = input_type
        metadata[:cond_type] = :none
        metadata[:product] = :dsgevarirf
        metadata[:class] = :obs # We default everything to an observable
        metadata[:date_inds] = OrderedDict()

        # Set up for loop over variable names
        means = DataFrame()
        bands = Dict{Symbol,DataFrame}()
        metadata[:indices] =
            OrderedDict{Symbol,Int}(name => name_i
                                    for (name_i, name) in enumerate(observables))

        # Means and Bands for each variable in a class
        for (name_i,name) in enumerate(observables)
            # irf_output is Vector{nperiod x nobs} -> for each observable,
            # we want to select its specific IRF, i.e. map(x -> x[:,obs_index]).
            # This creates a nperiod x ndraws matrix, which we want to transpose
            # to get a ndraws x nperiod matrix
            single_var = Matrix(reduce(hcat, map(x -> x[:,name_i], irf_output))')
            if correct_integ_series && name in [:laborproductivity, :Wages]
                zt_tmp = Matrix(reduce(hcat, map(x -> x[:,findfirst(:z_t .== observables)],
                                                 irf_output))')
                zt_tmp = cumsum(zt_tmp, dims = 2)
                single_var += zt_tmp
            end
            if annualize_inflation && name in [:obs_corepce, :obs_gdpdeflator, :π_t, :NominalWageGrowth, :Epi_t, :Eπ_t, :obs_longinflation]
                single_var = quartertoannual(single_var) # multiplies by 4
            end

            means[!,name] = vec(mean(single_var, dims = 1))
            bands[name]   = find_density_bands(single_var, density_bands;
                                               minimize = minimize)
        end

        if compute_inflation_expectation_error
            if :obs_gdpdeflator in observables
                gdpdef_i = findfirst(:obs_gdpdeflator .== observables)
                gdpdef_var = Matrix(reduce(hcat, map(x -> x[:, gdpdef_i], irf_output))')
                if :Epi_t in observables
                    espi_i = findfirst(:Epi_t .== observables)
                    espi_var = Matrix(reduce(hcat, map(x -> x[:, espi_i], irf_output))')
                    means[!,:Epi_err] = vec(mean(gdpdef_var - espi_var, dims = 1))
                    bands[:Epi_err]   = find_density_bands(gdpdef_var - espi_var,
                                                           density_bands;
                                                           minimize = minimize)
                elseif :obs_longinflation in observables
                    elpi_i = findfirst(:obs_longinflation .== observables)
                    elpi_var = Matrix(reduce(hcat, map(x -> x[:, elpi_i], irf_output))')
                    means[!,:long_Epi_err] = vec(mean(gdpdef_var - elpi_var, dims = 1))
                    bands[:long_Epi_err]   = find_density_bands(gdpdef_var - elpi_var,
                                                                density_bands;
                                                                minimize = minimize)
                end
            elseif :π_t in observables
                pi_i = findfirst(:π_t .== observables)
                pi_var = Matrix(reduce(hcat, map(x -> x[:, pi_i], irf_output))')
                if :Epi_t in observables
                    espi_i = findfirst(:Epi_t .== observables)
                    espi_var = Matrix(reduce(hcat, map(x -> x[:, espi_i], irf_output))')
                    means[!,:Epi_err] = vec(mean(pi_var - espi_var, dims = 1))
                    bands[:Epi_err]   = find_density_bands(pi_var - espi_var,
                                                           density_bands;
                                                           minimize = minimize)
                elseif :obs_longinflation in observables
                    elpi_i = findfirst(:obs_longinflation .== observables)
                    elpi_var = Matrix(reduce(hcat, map(x -> x[:, elpi_i], irf_output))')
                    means[!,:long_Epi_err] = vec(mean(pi_var - elpi_var, dims = 1))
                    bands[:long_Epi_err]   = find_density_bands(pi_var - elpi_var,
                                                                density_bands;
                                                                minimize = minimize)
                end
            end
        end
        mb = MeansBands(metadata, means, bands)

        # Save MeansBands
        tail = if method == :cholesky
            :_cholesky
        elseif method == :maxBC || method == :maximum_business_cycle_variance
            :_maxBC
        else
            :_choleskyLR
        end
        tail = Symbol(tail, addl_text)

        var_names = Symbol(join(string.(map(x -> DSGE.detexify(x), observables)), "_"))
        @show var_names
        fp = get_meansbands_output_file(m, input_type, :none,
                                        Symbol(:dsgevarirf, :obs_,
                                               var_names, tail),
                                        forecast_string = forecast_string)
        dirpath = dirname(fp)
        isdir(dirpath) || mkpath(dirpath)
        JLD2.jldopen(fp, true, true, true, IOStream) do file
            write(file, "mb", mb)
        end
        println(verbose, :high, "  " * "wrote " * basename(fp))
        return mb
    else
        # Reshape irf_output to nobs x nperiod x ndraw
        results = cat(map(x -> x', irf_output)..., dims = 3)
        if input_type == :mode
            for (k,var) in enumerate(observables)
                if annualize_inflation && var in [:obs_corepce, :obs_gdpdeflator, :π_t, :NominalWageGrowth, :Epi_t, :obs_longinflation]
                    results[k,:,:] = quartertoannual(results[k,:,:])
                end
            end
        end

        return results
    end
end
