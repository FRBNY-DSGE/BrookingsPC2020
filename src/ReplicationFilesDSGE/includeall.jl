# This file includes any scripts written specifically for the
# Brookings project on the Phillips Curve
using DataFrames, Dates, DSGE, DSGEModels, JLD2, KernelDensity, LinearAlgebra, MAT, ModelConstructors, Nullables, OrderedCollections
using Plots, Printf, Nullables, StatsPlots
using Plots.PlotMeasures, ColorTypes
using LinearAlgebra

include("util.jl")

# plot/
#include("plot/cond_obs_shocks_comparison_plot_forecast.jl")
#include("plot/plot_decomp_barplot_cond_obs_shocks.jl")
include("plot/plot_brookings_irfs.jl")
#include("plot/plot_shockdecirfbar_cond_obs_shocks.jl")
#include("plot/plot_shockdecirfline_cond_obs_shocks.jl")
#include("plot/shock_groups.jl")
include("plot/plot_4line.jl")

# forecast/
#include("forecast/decomp_cond_obs_shocks.jl")
include("forecast/cholesky_responses.jl")
include("forecast/impulse_responses.jl")
include("forecast/misc.jl")
include("forecast/util.jl")
include("forecast/dsgevar.jl")
include("forecast/dsgevar_impulse_responses.jl")
#include("forecast/dsgevar_estimate.jl")
include("forecast/var_impulse_responses.jl")

# scenarios/
#include("scenarios/brookingsPC_scenarios.jl")

# write_tex/
#include("write_tex/write_cond_obs_shocks_comparison.jl")
#include("write_tex/write_brookings_irf.jl")
#include("write_tex/util.jl")
