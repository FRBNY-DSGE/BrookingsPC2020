using DSGE, OrderedCollections
using ClusterManagers, HDF5
import DSGE: usual_settings!, usual_forecast!
using ModelConstructors, Dates

include("../includeall.jl")

estimate = true
sample = "incZLB"
test_loading_in = false

# Initialize model objects
m = Model1002("ss21") #, custom_settings = custom_settings)
fp = dirname(@__FILE__)
standard_spec!(m, "191118", fp, fcast_date = Date(2019, 12, 31))
m <= Setting(:add_laborshare_measurement, false)
m <= Setting(:hours_first_observable, false)
m <= Setting(:data_vintage, "191118")
m <= Setting(:date_forecast_start, quartertodate("2019-Q4"))


m <= Setting(:saveroot, joinpath(fp, "../../../save/"))
m <= Setting(:dataroot, joinpath(fp, "../../../save/input_data"))

m <= Setting(:date_regime2_start, Date(1990, 3, 31))
m <= Setting(:date_regime2_start_text, "900331", true, "reg2start", "The text version to be saved of when regime 2 starts")

df = load_data(m)

m <= Setting(:n_particles, 15000)
m <= Setting(:n_smc_blocks, 5)
m <= Setting(:n_mh_steps_smc, 1)
m <= Setting(:use_parallel_workers, true)
m <= Setting(:step_size_smc, 0.5)
m <= Setting(:resampler_smc, :systematic)
m <= Setting(:target_accept, 0.25)
m <= Setting(:mixture_proportion, 0.9)
m <= Setting(:use_fixed_schedule, false)
m <= Setting(:adaptive_tempering_target_smc, 0.97)
m <= Setting(:resampling_threshold, .5)
m <= Setting(:use_fixed_schedule, false)

if estimate
    df_full_preZLB = df[Date(1964, 3, 1) .<= df[:date] .<= get_setting(m, :date_zlb_start) - Dates.Month(3), :]
    df_full_incZLB = df[Date(1964, 3, 1) .<= df[:date], :]

  #  my_procs = addprocs_frbny(100)
    @everywhere using DSGE, OrderedCollections

    if sample=="preZLB"
        m <= Setting(:period, "full_preZLB", true, "period", "period for this exercse (first or second)")
        m <= Setting(:n_regimes, 2, true, "reg", "Number of Regimes")
        m <= Setting(:date_presample_start, DSGE.quartertodate("1964-Q1")) # need to reset presample
        m <= Setting(:date_mainsample_start, DSGE.quartertodate("1964-Q3"))
        smc2(m, df_to_matrix(m, df_full_preZLB), save_intermediate = true, intermediate_stage_increment = 1, regime_switching = true, n_regimes = 2)
    elseif sample=="incZLB"
        m <= Setting(:period, "full_incZLB", true, "period", "period for this exercse (first or second)")
        m <= Setting(:n_regimes, 2, true, "reg", "Number of Regimes")
        m <= Setting(:date_presample_start, DSGE.quartertodate("1964-Q1")) # need to reset presample
        m <= Setting(:date_mainsample_start, DSGE.quartertodate("1964-Q3"))
        smc2(m, df_to_matrix(m, df_full_incZLB), save_intermediate = true, intermediate_stage_increment = 1, regime_switching = true, n_regimes = 2)
    end

    rmprocs(my_procs)

end

if test_loading_in
    using JLD2, FileIO, Statistics
    cloud = load("../../../save/output_data/$(spec(m))/$(subspec(m))/estimate/raw/smc_cloud_period=full_preZLB_reg2start=900331_reg=2_vint=191118.jld2", "cloud")
    draws = DSGE.get_vals(cloud)
    @show mean(draws[2,:])
    @show mean(draws[21, :])
    @show quantile(draws[2,:], [0.05, 0.95])
    @show quantile(draws[21,:], [0.05, 0.95])
    aaa
    overrides = forecast_input_file_overrides(m)
    overrides[:full] = joinpath(fp, "../../../save/output_data/$(spec(m))/$(subspec(m))/estimate/raw/smcsave_period=full_preZLB_reg=2_vint=191118.h5")
    draws = load_draws(m, :full)
end
