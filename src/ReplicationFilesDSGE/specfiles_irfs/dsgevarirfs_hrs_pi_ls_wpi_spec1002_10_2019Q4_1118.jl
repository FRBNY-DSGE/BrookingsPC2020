using DSGE, ModelConstructors, Dates, LinearAlgebra, Test
include("../includeall.jl")
# ENV["frbnyjuliamemory"] = "1G"

## What to do?
run_irfs         = false
do_rev_transform = false
do_add_workers   = false
make_plots       = true
n_workers = 50

## Set flags
vint = "191118"
forecast_string1 = "first"
forecast_string2 = "second"
obs_varnames    = Vector{Symbol}(undef,0)
pseudo_varnames = Vector{Symbol}(undef,0)
push!(obs_varnames, [:obs_hours, :obs_gdp, :obs_gdpdeflator, :obs_wages]...)
push!(pseudo_varnames, [:Hours, :y_t, :π_t, :Wages, :laborshare_t, :labor_productivity]...)
density_bands = [.9]

# Set up DSGE-VAR system
lags = 4
# observables = [:LaborProductivityGrowthNoME, :Hours, :obs_corepce,
#                :laborshare_t, :NominalWageGrowth]
# shocks = [:zp_sh, :b_sh, :rm_sh, :λ_f_sh, :ztil_sh]
shocks = Vector{Symbol}(undef,0) # empty implies all shocks
observables = [:obs_hours, :π_t, :laborshare_t, :NominalWageGrowth]
obs_shock = :obs_hours

# Set up ylims for plots
ylims_dict = Dict{Symbol,Tuple{Float64,Float64}}()
ylims_dict[:obs_hours] = (-2.75, -0.5)
ylims_dict[:π_t] = (-0.6, 0.2)
ylims_dict[:laborshare_t] = (-1.0, 0.2)
ylims_dict[:NominalWageGrowth] = (-1.3, 0.3)

ylabels_dict = Dict{Symbol,String}()
ylabels_dict[:obs_hours] = "percent"
ylabels_dict[:π_t] = "percentage points"
ylabels_dict[:laborshare_t] = "percent"
ylabels_dict[:NominalWageGrowth] = "percentage points"

xlabels_dict = Dict{Symbol,String}()
xlabels_dict[:obs_hours] = "horizon"
xlabels_dict[:π_t] = "horizon"
xlabels_dict[:laborshare_t] = "horizon"
xlabels_dict[:NominalWageGrowth] = "horizon"

# ylims_dict[:laborproductivity] = (-0.4, 0.5)
# ylims_dict[:MarginalCost] = (-0.1, 1.0)

## Set up model and data for subsamples
# dsid is "model spec"sample number, e.g. 10021
# cdid indicates which variables ordering for
# conditional forecast.
# cdid = 1 => hours first and 1pp increase
custom_settings1 = Dict{Symbol,Setting}(:add_laborshare_measurement =>
                                        Setting(:add_laborshare_measurement, true),
                                        :hours_first_observable =>
                                        Setting(:hours_first_observable, false),
                                        :add_laborproductivity_measurement =>
                                        Setting(:add_laborproductivity_measurement, false),
                                        :add_laborproductivitygrowth_nome_measurement =>
                                        Setting(:add_laborproductivitygrowth_nome_measurement,
                                                false))
custom_settings2 = Dict{Symbol,Setting}(:add_laborshare_measurement =>
                                        Setting(:add_laborshare_measurement, true),
                                        :hours_first_observable =>
                                        Setting(:hours_first_observable, false),
                                        :add_laborproductivity_measurement =>
                                        Setting(:add_laborproductivity_measurement, false),
                                        :add_laborproductivitygrowth_nome_measurement =>
                                        Setting(:add_laborproductivitygrowth_nome_measurement,
                                                false))

fp = dirname(@__FILE__)
m1 = Model1002("ss10"; custom_settings = custom_settings1)
standard_spec200129!(m1, vint, fp; fcast_date = Date(2019, 12, 31), dsid = 10021, cdid = 1)
m1 <= Setting(:period, "r1", true, "period", "period for this exercse (first or second)")
m1 <= Setting(:preZLB, "false", true, "preZLB", "")
m1 <= Setting(:npart, "15000", true, "npart", "")
m1 <= Setting(:use_population_forecast, false)
m2 = Model1002("ss10"; custom_settings = custom_settings2)
standard_spec200129!(m2, vint, fp; fcast_date = Date(2019, 12, 31), dsid = 10022, cdid = 1)
m2 <= Setting(:period, "r2", true, "period", "period for this exercse (first or second)")
m2 <= Setting(:preZLB, "false", true, "preZLB", "")
m2 <= Setting(:npart, "15000", true, "npart", "")
m2 <= Setting(:use_population_forecast, false)
m1 <= Setting(:date_regime2_start_text, "900331", true, "reg2start",
              "The text version to be saved of when regime 2 starts")
m2 <= Setting(:date_regime2_start_text, "900331", true, "reg2start",
              "The text version to be saved of when regime 2 starts")

try
    @assert !isempty(shocks)
catch
    push!(shocks, collect(keys(m1.exogenous_shocks))...)
end

if run_irfs
    if do_add_workers
        addprocs_frbny(n_workers)
        @everywhere using DSGE, DSGEModels, OrderedCollections
    end

    parallel = do_add_workers ? true : false
    # impulse_responses(m1, load_draws(m1, :mode), :mode, :cholesky,
    #                   lags, observables, shocks,
    #                   findfirst(observables .== obs_shock);
    #                   parallel = parallel, flip_shocks = false,
    #                   forecast_string = forecast_string1)
    # impulse_responses(m2, load_draws(m2, :mode), :mode, :cholesky,
    #                   lags, observables, shocks,
    #                   findfirst(observables .== obs_shock);
    #                   parallel = parallel, flip_shocks = false,
    #                   forecast_string = forecast_string2)
    impulse_responses(m1, load_draws(m1, :full), :full, :cholesky,
                      lags, observables, shocks,
                      findfirst(observables .== obs_shock);
                      parallel = parallel, flip_shocks = false,
                      compute_meansbands = true,
                      forecast_string = forecast_string1)
    impulse_responses(m2, load_draws(m2, :full), :full, :cholesky,
                      lags, observables, shocks,
                      findfirst(observables .== obs_shock);
                      parallel = parallel, flip_shocks = false,
                      compute_meansbands = true,
                      forecast_string = forecast_string2)
    impulse_responses(m1, load_draws(m1, :full), :full, :maxBC,
                      lags, observables, shocks,
                      findfirst(observables .== obs_shock);
                      parallel = parallel, flip_shocks = false,
                      compute_meansbands = true,
                      forecast_string = forecast_string1)
    impulse_responses(m2, load_draws(m2, :full), :full, :maxBC,
                      lags, observables, shocks,
                      findfirst(observables .== obs_shock);
                      parallel = parallel, flip_shocks = false,
                      compute_meansbands = true,
                      forecast_string = forecast_string2)

    @info "Completed calculation of IRFs"
end

# Plot IRFs
if make_plots
    obs_ind    = findall([k in collect(keys(m1.observables)) for k in observables])
    pseudo_ind = findall([k in collect(keys(m1.pseudo_observables)) for k in observables])
    outobs = plot_brookings_irfs(m1, m2, observables[obs_ind],
                                 :obs, :cholesky, :full, :none;
                                 irf_product = :dsgevarirf, shocks = shocks, lags = lags, obs_shock = obs_shock,
                                 observables = observables,
                                 forecast_string1 = forecast_string1,
                                 forecast_string2 = forecast_string2,
                                 legend_names = ["Pre 1990", "Post 1990"],
                                 ylims_dict = ylims_dict,
                                 ylabels_dict = ylabels_dict,
                                 xlabels_dict = xlabels_dict)
    outpse = plot_brookings_irfs(m1, m2, observables[pseudo_ind],
                        :pseudo, :cholesky, :full, :none;
                        irf_product = :dsgevarirf, shocks = shocks, lags = lags, obs_shock = obs_shock,
                        observables = observables,
                        forecast_string1 = forecast_string1,
                        forecast_string2 = forecast_string2,
                        legend_names = ["Pre 1990", "Post 1990"],
                        ylims_dict = ylims_dict,
                                 ylabels_dict = ylabels_dict,
                                 xlabels_dict = xlabels_dict)
    plot_brookings_irfs(m1, m2, observables[obs_ind],
                        :obs, :maxBC, :full, :none;
                        irf_product = :dsgevarirf, shocks = shocks, lags = lags, obs_shock = obs_shock,
                        observables = observables,
                        forecast_string1 = forecast_string1,
                        forecast_string2 = forecast_string2,
                        legend_names = ["Pre 1990", "Post 1990"],
                        ylims_dict = ylims_dict,
                                 ylabels_dict = ylabels_dict,
                                 xlabels_dict = xlabels_dict)
    plot_brookings_irfs(m1, m2, observables[pseudo_ind],
                        :pseudo, :maxBC, :full, :none;
                        irf_product = :dsgevarirf, shocks = shocks, lags = lags, obs_shock = obs_shock,
                        observables = observables,
                        forecast_string1 = forecast_string1,
                        forecast_string2 = forecast_string2,
                        legend_names = ["Pre 1990", "Post 1990"],
                        ylims_dict = ylims_dict,
                                 ylabels_dict = ylabels_dict,
                                 xlabels_dict = xlabels_dict)

    # plot_brookings_irfs(m1, m2, observables[obs_ind],
    #                     :obs, :choleskyLR, :full, :none;
    #                     irf_product = :dsgevarirf, shocks = shocks, lags = lags, obs_shock = obs_shock,
    #                     observables = observables,
    #                     forecast_string1 = forecast_string1,
    #                     forecast_string2 = forecast_string2,
    #                     legend_names = ["Pre 1990", "Post 1990"],
    #                     ylims_dict = ylims_dict)
    # plot_brookings_irfs(m1, m2, observables[pseudo_ind],
    #                     :pseudo, :choleskyLR, :full, :none;
    #                     irf_product = :dsgevarirf, shocks = shocks, lags = lags, obs_shock = obs_shock,
    #                     observables = observables,
    #                     forecast_string1 = forecast_string1,
    #                     forecast_string2 = forecast_string2,
    #                     legend_names = ["Pre 1990", "Post 1990"],
    #                     ylims_dict = ylims_dict)

    # # Create shock decompositions at mode
    # DSGE.update!(m1, load_draws(m1, :mode))
    # DSGE.update!(m2, load_draws(m2, :mode))
    # system1 = compute_system(m1)
    # system2 = compute_system(m2)
    # permute_mat = Matrix{Float64}(I, n_observables(m1), n_observables(m1))
    # h = impulse_response_horizons(m1)
    # _, obs1, _ = impulse_responses(m1, system1, h, permute_mat,
    #                                vcat(1.,
    #                                     zeros(n_observables(m1) - 1)),
    #                                flip_shocks = false,
    #                                cholesky_obs_shock = true)
    # hr_shock1 = 1. / obs1[1,1]
    # states1, obs1, pseudo1, _, irfshock1 = impulse_responses(m1, system1, h, permute_mat,
    #                                           vcat(hr_shock1,
    #                                                zeros(n_observables(m1) - 1)),
    #                                           cholesky_obs_shock = true,
    #                                           flip_shocks = false,
    #                                           get_shocks = true)
    # _, obs2, _ = impulse_responses(m2, system2, h, permute_mat,
    #                                vcat(1.,
    #                                     zeros(n_observables(m2) - 1)),
    #                                flip_shocks = false,
    #                                cholesky_obs_shock = true)
    # hr_shock2 = 1. / obs2[1,1]
    # _, obs2, pseudo2, _, irfshock2 = impulse_responses(m2, system2, h, permute_mat,
    #                                           vcat(hr_shock2,
    #                                                zeros(n_observables(m2) - 1)),
    #                                           flip_shocks = false,
    #                                           cholesky_obs_shock = true,
    #                                           get_shocks = true)

    # if do_rev_transform
    #     for (k,v) in m1.observable_mappings
    #         irf_trans = DSGE.get_irf_transform(v.rev_transform)
    #         obs1[m1.observables[k],:] = irf_trans(obs1[m1.observables[k],:])
    #     end
    #     for (k,v) in m2.observable_mappings
    #         irf_trans = DSGE.get_irf_transform(v.rev_transform)
    #         obs2[m2.observables[k],:] = irf_trans(obs2[m2.observables[k],:])
    #     end
    # end

    # pseudo1 = correct_integ_series!(m1, pseudo1)
    # pseudo2 = correct_integ_series!(m2, pseudo2)
    # df_irfshock1 = DataFrame()
    # df_irfshock2 = DataFrame()
    # df_obs1 = DataFrame()
    # df_obs2 = DataFrame()
    # df_pseudo1 = DataFrame()
    # df_pseudo2 = DataFrame()
    # for (shock,ind) in m1.exogenous_shocks
    #     df_irfshock1[!,shock] = [irfshock1[ind]]
    #     df_irfshock2[!,shock] = [irfshock2[ind]]
    # end
    # for (k,v) in m1.observables
    #     df_obs1[!,k] = vcat(0., obs1[v,:])
    #     df_obs2[!,k] = vcat(0., obs2[v,:])
    # end
    # for (k,v) in m1.pseudo_observables
    #     df_pseudo1[!,k] = vcat(0., pseudo1[v,:])
    #     df_pseudo2[!,k] = vcat(0., pseudo2[v,:])
    # end
    # plot_shockdecirfline_cond_obs_shocks(m1, :mode, :none, [:pseudo, :obs],
    #                                      shock_groupings(m1),
    #                                      df_irfshock1,
    #                                      [DataFrame(df_pseudo1[2:end,:]),
    #                                       DataFrame(df_obs1[2:end,:])],
    #                                      forecast_string = forecast_string1,
    #                                      forecast_string1 = forecast_string1,
    #                                      forecast_string2 = forecast_string1,
    #                                      fourquarter = [false, false],
    #                                      addtl_fn_text = "_cholesky",
    #                                      ylims_dict = ylims_dict)
    # plot_shockdecirfline_cond_obs_shocks(m2, :mode, :none, [:pseudo, :obs],
    #                                      shock_groupings(m2), df_irfshock2,
    #                                      [DataFrame(df_pseudo2[2:end,:]),
    #                                       DataFrame(df_obs2[2:end,:])],
    #                                      forecast_string = forecast_string2,
    #                                      forecast_string1 = forecast_string2,
    #                                      forecast_string2 = forecast_string2,
    #                                      fourquarter = [false, false],
    #                                      addtl_fn_text = "_cholesky",
    #                                      ylims_dict = ylims_dict)
end

nothing
