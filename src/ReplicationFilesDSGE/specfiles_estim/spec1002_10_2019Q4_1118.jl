using DSGE, OrderedCollections
using ClusterManagers, HDF5
import DSGE: usual_settings!, usual_forecast!
using ModelConstructors, Dates

# Initialize model objects
m = Model1002("ss10")

m <= Setting(:data_vintage, "191118")
fp = dirname(@__FILE__)
m <= Setting(:saveroot, joinpath(fp, "../../../save/"))
m <= Setting(:dataroot, joinpath(fp, "../../../save/input_data/"))
m <= Setting(:date_forecast_start, Date(2019,12,31))

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

df_first = df[Date(1964, 3, 1) .<= df[:date]  .<= Date(1989, 12, 31), :]
df_second = df[Date(1990, 3, 1) .<= df[:date], :]


my_procs = addprocs(200)
@everywhere using DSGE, OrderedCollections, DSGEModels

# FIRST
m <= Setting(:period, "first", true, "period", "period for this exercse (first or second)")
m <= Setting(:date_presample_start, DSGE.quartertodate("1964-Q1")) # need to reset presample
m <= Setting(:date_mainsample_start, DSGE.quartertodate("1964-Q3"))
smc2(m, df_to_matrix(m, df_first))

# SECOND
m <= Setting(:period, "second", true, "period", "period for this exercse (first or second)")
m <= Setting(:date_presample_start, DSGE.quartertodate("1990-Q1")) # need to reset presample
m <= Setting(:date_mainsample_start, DSGE.quartertodate("1990-Q3"))
smc2(m, df_to_matrix(m, df_second))

rmprocs(my_procs)
