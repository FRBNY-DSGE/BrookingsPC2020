using DSGE, ModelConstructors, Dates, LinearAlgebra, Test, NLsolve, Roots
include("../includeall.jl")


save_orig = true

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
models = Dict{Symbol,AbstractDSGEModel}()
for model in [:H1, :κp, :MPparam]
    m = Model1002("ss10"; custom_settings = custom_settings1)
    standard_spec!(m, vint, fp; fcast_date = Date(2019, 12, 31), dsid = 10021, cdid = 1, save_orig = save_orig)
    m <= Setting(:period, "r1", true, "period", "period for this exercse (first or second)")
    m <= Setting(:npart, "15000", true, "npart", "")
    m <= Setting(:preZLB, false, true, "preZLB", "")
    m <= Setting(:reg2start, "900331", true, "reg2start", "")
    m <= Setting(:use_population_forecast, false)
    DSGE.update!(m, load_draws(m, :mode; use_highest_posterior_value = true))
    models[model] = m
end
showparams = [:ζ_p, :ψ1, :ψ2, :ψ3, :ρ, :ι_p, :ι_w]

# Get H2
m2 = Model1002("ss10"; custom_settings = custom_settings2)
standard_spec!(m2, vint, fp; fcast_date = Date(2019, 12, 31), dsid = 10022, cdid = 1, save_orig = save_orig)
m2 <= Setting(:period, "r2", true, "period", "period for this exercse (first or second)")
m2 <= Setting(:npart, "15000", true, "npart", "period for this exercse (first or second)")
m2 <= Setting(:preZLB, false, true, "preZLB", "")
m2 <= Setting(:reg2start, "900331", true, "reg2start", "")
m2 <= Setting(:use_population_forecast, false)
DSGE.update!(m2, load_draws(m2, :mode; use_highest_posterior_value = true))
models[:H2] = m2

# Update mp rule coefficients
θ_MPparam = map(x -> x.value, models[:MPparam].parameters)
θ_MPparam[models[:MPparam].keys[:ρ]] = m2[:ρ].value
θ_MPparam[models[:MPparam].keys[:ψ1]] = m2[:ψ1].value
θ_MPparam[models[:MPparam].keys[:ψ2]] = m2[:ψ2].value
θ_MPparam[models[:MPparam].keys[:ψ3]] = m2[:ψ3].value
DSGE.update!(models[:MPparam], θ_MPparam)

# Show H1 vs H2 parameters
for i in showparams
    println("Showing parameter: " * string(i))
    println("($(models[:H1][i].value), $(models[:H2][i].value))")
end
κp1 = get_κp(models[:H1], reshape(map(x -> x.value, models[:H1].parameters),
                                  1, n_parameters(models[:H1])))
κw1 = get_κw(models[:H1], reshape(map(x -> x.value, models[:H1].parameters),
                                  1, n_parameters(models[:H1])))

κp2 = get_κp(models[:H2], reshape(map(x -> x.value, models[:H2].parameters),
                                  1, n_parameters(models[:H2])))
κw2 = get_κw(models[:H2], reshape(map(x -> x.value, models[:H2].parameters),
                                  1, n_parameters(models[:H2])))
println("Showing parameter: κp")
println("($(κp1), $(κp2))")
println("Showing parameter: κw")
println("($(κw1), $(κw2))")

# Find the ζp which makes κp1 = κp2
function f!(F, x)
    ζ_p = x
    θ = map(x -> x.value, models[:κp].parameters)
    θ[models[:κp].keys[:ζ_p]] = ζ_p[1]
    DSGE.update!(models[:κp], θ)
    F = κp2 - get_κp(models[:κp], reshape(map(x -> x.value, models[:κp].parameters), 1,
                                          n_parameters(models[:κp])))
end

function f(x)
    ζ_p = x
    θ = map(x -> x.value, models[:κp].parameters)
    θ[models[:κp].keys[:ζ_p]] = ζ_p
    DSGE.update!(models[:κp], θ)
    return κp2 - get_κp(models[:κp], reshape(map(x -> x.value, models[:κp].parameters), 1,
                                          n_parameters(models[:κp])))
end

models[:κp] <= Setting(:fix_ζ_p, models[:κp][:ζ_p].value)
out2 = find_zero(f, (1e-6, 1-1e-6), Bisection())
θ = map(x -> x.value, models[:κp].parameters)
θ[models[:κp].keys[:ζ_p]] = out2
DSGE.update!(models[:κp], θ)
@show out2
@show get_κp(models[:κp], reshape(map(x -> x.value, models[:κp].parameters), 1,
                                n_parameters(models[:κp])))
@show κp2
try
    @assert !isempty(shocks)
catch
    push!(shocks, collect(keys(m2.exogenous_shocks))...)
end

cholesky_irf = Dict()
maxBC_irf = Dict()
for (k,mod) in models
    cholesky_irf[k] =
        impulse_responses(mod, reshape(map(x -> x.value, mod.parameters),
                                     1, length(mod.parameters)),
                          :mode, :cholesky,
                          lags, observables, shocks,
                          findfirst(observables .== obs_shock);
                          parallel = false, flip_shocks = false)
    maxBC_irf[k] =
        impulse_responses(mod, reshape(map(x -> x.value, mod.parameters),
                                     1, length(mod.parameters)),
                          :mode, :maxBC,
                          lags, observables, shocks,
                          findfirst(observables .== obs_shock);
                          parallel = false, flip_shocks = false)
end

@info "Completed calculation of IRFs"

# Plot IRFs
df_obss_cholesky = OrderedDict{Symbol,DataFrame}()
df_obss_maxBC = OrderedDict{Symbol,DataFrame}()
for k in keys(models)
    df_obss_cholesky[k] = DataFrame()
    df_obss_maxBC[k] = DataFrame()
end
for v in observables
    for (key, m) in models
        df_obss_cholesky[key][!,v] = vec(cholesky_irf[key][findfirst(v .== observables),:,:])
        df_obss_maxBC[key][!,v] = vec(maxBC_irf[key][findfirst(v .== observables),:,:])
    end
end

plot_dict = Dict{Symbol,Dict{Symbol,Plots.Plot}}()
for (dfs, pname) in zip([df_obss_cholesky, df_obss_maxBC], [:cholesky, :maxBC])
    plot_dict[pname] = Dict{Symbol, Plots.Plot}()
    for var in observables
        detexify_title = typeof(Plots.backend()) == Plots.GRBackend
        if var in collect(keys(models[:H1].observables))
            the_title = DSGE.describe_series(models[:H1], var, :obs, detexify = detexify_title)
        else
            the_title = DSGE.describe_series(models[:H1], var, :pseudo, detexify = detexify_title)
        end
        p = plot(1:20, dfs[:H1][!,var], label = "Pre 1990", color = RGB(55. / 255., 126. / 255., 184. / 255.),
                 legend = :bottomright, linewidth = 2,
                 xtickfont = font(12), ytickfont = font(12),
                 left_margin = 50px, right_margin = 10px,
                 bottom_margin = 40px, xlabel = "horizon", ylabel = ylabels_dict[var])
        plot!(p, 1:20, dfs[:H2][!,var], label = "Post 1990", color = RGB(.8941, .1020, .1098),
              linewidth = 2)
        plot!(p, 1:20, dfs[:κp][!,var], label = "Counterfactual Slope", color = :black,
              linewidth = 2)
        plot!(p, 1:20, dfs[:MPparam][!,var], label = "Counterfactual Policy", color = :black,
              linestyle = :dash,
              linewidth = 2, ylims = ylims_dict[var])
        plot_dict[pname][var] = p

        savefig(p, figurespath(models[:H1], "forecast", "counter_fix_kappap_dsgevarirf_" * string(DSGE.detexify(var)) * "_" * string(pname) * "_variables_hrs_pi_ls_wpi.pdf"))
    end
end
nothing
