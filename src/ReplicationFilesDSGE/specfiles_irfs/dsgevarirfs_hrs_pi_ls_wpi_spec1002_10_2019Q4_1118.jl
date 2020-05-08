using DSGE, ModelConstructors, Dates, LinearAlgebra, Test

save_orig = true

include("../includeall.jl")

## What to do?
run_irfs         = true
do_rev_transform = false
do_add_workers   = true
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

## Set up model and data for subsamples
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
standard_spec!(m1, vint, fp; fcast_date = Date(2019, 12, 31), dsid = 10021, cdid = 1, save_orig = save_orig)
m1 <= Setting(:period, "r1", true, "period", "period for this exercse (first or second)")
m1 <= Setting(:preZLB, "false", true, "preZLB", "")
m1 <= Setting(:npart, "15000", true, "npart", "")
m1 <= Setting(:use_population_forecast, false)
m2 = Model1002("ss10"; custom_settings = custom_settings2)
standard_spec!(m2, vint, fp; fcast_date = Date(2019, 12, 31), dsid = 10022, cdid = 1, save_orig = save_orig)
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
        @everywhere include("../includeall.jl")
    end

    parallel = do_add_workers ? true : false
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
end

nothing
