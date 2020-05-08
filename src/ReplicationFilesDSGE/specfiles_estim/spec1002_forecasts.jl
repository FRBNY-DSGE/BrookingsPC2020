using DSGE, Dates, ModelConstructors, Plots, Plots.PlotMeasures, FileIO
using Nullables
include("../util.jl")
include("../ait.jl")
include("../plt.jl")
include("../zero_pi_star.jl")

mode = "smc"
gap_value = 0.0

m = Model1002("ss21")
est_vint = "191118"
data_vint = "200227"

standard_spec!(m, data_vint)
m <= Setting(:data_vintage, data_vint)
m <= Setting(:date_forecast_start, quartertodate("2020-Q1"))

m <= Setting(:saveroot, "../../../save_orig/")
m <= Setting(:dataroot, "../../../save/input_data")

m <= Setting(:date_regime2_start, Date(1990, 3, 31))
m <= Setting(:date_regime2_start_text, "900331", true, "reg2start", "The text version to be saved of when regime 2 starts")

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

m <= Setting(:period, "full_incZLB", true, "period", "period for this exercse (first or second)")
m <= Setting(:n_regimes, 2, true, "reg", "Number of Regimes")
m <= Setting(:regime_switching, true)
m <= Setting(:date_presample_start, DSGE.quartertodate("1964-Q1")) # need to reset presample
m <= Setting(:date_mainsample_start, DSGE.quartertodate("1964-Q3"))

output_vars = [:forecastobs, :histobs, :forecastpseudo, :histpseudo]
data = load_data(m)
param_symbols = map(x->x.key, m.parameters)
ind = findfirst(x -> x==:ζ_p, param_symbols)
ind2 = findfirst(x -> x==:ζ_p_r2, param_symbols)
ind_psi_pi = findfirst(x -> x==:ψ1, param_symbols)

m <= Setting(:data_vintage, est_vint)
param_mode = load_draws(m, :mode)

savepath = saveroot(m)

mode_use = if mode == "csminwel"
    if !@isdefined csminwel_mode
        csminwel_mode = load("csminwel.jld2", "csminwel_mode")
    end
    csminwel_mode
elseif mode == "smc"
    param_mode
end

plot_start_date = Date(2007, 3, 31)

function new_forecast_model(yr::Int)
    m_forecast = Model1002("ss13")
    m_forecast <= Setting(:date_presample_start, DSGE.quartertodate("$(string(yr))-Q1"))
    m_forecast <= Setting(:date_mainsample_start, DSGE.quartertodate("$(string(yr))-Q3"))
    m_forecast <= Setting(:forecast_horizons, 32)
    return m_forecast
end


for yr in [1990]
    for obs in [:obs_corepce, :MarginalCost]
        for e in [""]
            plots = Dict{Symbol, Plots.Plot}()
            plots[:plt] = plot()
            plots[:ait] = plot()
            plots[:psi] = plot()

            m_forecast = new_forecast_model(yr)

            r2_params = vcat(mode_use[1:(ind-1)], mode_use[ind2],
                             mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end])
            r2_kappa_p = get_κp(m_forecast, r2_params)
            r1_params = mode_use[1:end .!=ind2]
            r1_kappa_p = get_κp(m_forecast, r1_params)
            ind_eps_p = findfirst(x -> x==:ϵ_p, map(x->x.key, m_forecast.parameters))
            while abs(r2_kappa_p - r1_kappa_p) >= 0.00001
                m_forecast[:ϵ_p] = m_forecast[:ϵ_p] + 0.001
                r1_kappa_p = get_κp(m_forecast, r1_params)
            end
            eps_val = m_forecast[:ϵ_p].value

            m_forecast = new_forecast_model(yr)

            # Forecast using R1 (Blue)
            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars, mode_use[1:end .!=ind2], data)
            if obs in [:MarginalCost, :Sinf_t, :πtil_t]
                last_val = fcast[:histpseudo][m_forecast.pseudo_observables[obs], :][end]

                qtrs = size(fcast[:histpseudo], 2) + size(fcast[:forecastpseudo], 2) #- 1
                start_qtr = Date(yr, 9, 30)
                dates = start_qtr:Dates.Month(3):iterate_quarters(start_qtr, qtrs)
                dates = map(x -> "$(year(x))-Q$(quarterofyear(x))", dates)

                mult = if obs == :πtil_t
                    4.0
                else
                    1.0
                end

                sys = compute_system(m_forecast)
                irf_s, irf_o, irf_p = impulse_responses(m_forecast, sys)
                irf_plot_b_sh = plot(irf_p[m_forecast.pseudo_observables[obs], :,
                                           m_forecast.exogenous_shocks[:b_sh]],
                                     title = "$(string(DSGE.detexify(obs))) b shock", label = "R1",
                                     left_margin = 30px, color = :blue)

                irf_plot_rm_sh = plot(irf_p[m_forecast.pseudo_observables[obs], :,
                                            m_forecast.exogenous_shocks[:rm_sh]],
                                      title = "$(string(DSGE.detexify(obs))) rm shock", label = "R1",
                                      left_margin = 30px, color = :blue)

            ### REGIME 2 PARAMETER VALUE (Red)
            m_forecast = Model1002("ss13")
            m_forecast <= Setting(:forecast_horizons, 32)
            m_forecast <= Setting(:alternative_policy, AltPolicy(:historical, eqcond, solve, forecast_init = identity))
            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data)

            last_val = fcast[:histpseudo][m_forecast.pseudo_observables[obs], :][end]

            for p in keys(plots)
                plot!(plots[p], dates[end-84:end],
                      mult*vcat(fcast[:histpseudo][m_forecast.pseudo_observables[obs], :], repeat([NaN], length(fcast[:forecastpseudo][m_forecast.pseudo_observables[obs], :])))[end-84:end],
                      label = "", color = :black)
                plot!(plots[p], dates[end-84:end], mult*vcat(repeat([NaN], length(fcast[:histpseudo][m_forecast.pseudo_observables[obs], :])-1), last_val, fcast[:forecastpseudo][m_forecast.pseudo_observables[obs], :])[end-84:end], label = "Post 1990 Slope", color = :red,
                      xticks = (1:20:length(dates[end-84:end]), dates[end-84:end][1:20:length(dates[end-84:end])]))
            end

            # REGIME 2 PARAMETER VALUE. RAISE PSI_PI
            m_forecast <= Setting(:alternative_policy, AltPolicy(:zero_pi_star, eqcond_zero_pi_star,
                                                                 solve_zero_pi_star,
                                                                 forecast_init = forecast_init_zero_pi_star))
            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data,
                                           param_key = :ψ1, param_value = 300.)
            last_val = fcast[:histpseudo][m_forecast.pseudo_observables[obs], :][end]
            plot!(plots[:psi], dates[end-84:end], mult*vcat(repeat([NaN], length(fcast[:histpseudo][m_forecast.pseudo_observables[obs], :])), fcast[:forecastpseudo][m_forecast.pseudo_observables[obs], :])[end-84:end], label = "R2, psi_pi = 300", color = :red, linestyle = :dash)

            m_forecast = Model1002("ss13")
            m_forecast <= Setting(:forecast_horizons, 32)
            sys = compute_system(m_forecast)
            irf_s, irf_o, irf_p = impulse_responses(m_forecast, sys)
            plot!(irf_plot_b_sh, irf_p[m_forecast.pseudo_observables[obs], :, m_forecast.exogenous_shocks[:b_sh]], title = "$(string(DSGE.detexify(obs))) b shock", label = "R2 zetap", color = :red)

            plot!(irf_plot_rm_sh, irf_p[m_forecast.pseudo_observables[obs], :, m_forecast.exogenous_shocks[:rm_sh]], title = "$(string(DSGE.detexify(obs))) rm shock", label = "R2 zetap", color = :red)


            ### REGIME 1 PARAMETER VALUE (Blue) #i use ind2 in the vcat since that's for the history, but ind in the param_value since that's for the forecast
            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data, param_key = :ζ_p, param_value = mode_use[ind])
            nans = repeat([NaN], length(fcast[:histpseudo][m_forecast.pseudo_observables[obs], :]))
            for p in keys(plots)
                plot!(plots[p], dates[end-84:end], mult*vcat(nans,
                                                   last_val,
                                                   fcast[:forecastpseudo][m_forecast.pseudo_observables[obs], :])[end-84:end], label = "Pre 1990 Slope", color = :blue)
            end


            ### REGIME 1 PARAMETER VALUE. Raise Ψ_π
            m_forecast <= Setting(:alternative_policy, AltPolicy(:zero_pi_star, eqcond_zero_pi_star,
                                                                 solve_zero_pi_star,
                                                                 forecast_init = forecast_init_zero_pi_star))
            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data, param_key = :ζ_p, param_value = mode_use[ind],
                                           param_key2 = :ψ1, param_value2 = 300.)
            nans = repeat([NaN], length(fcast[:histpseudo][m_forecast.pseudo_observables[obs], :]))
            plot!(plots[:psi], dates[end-84:end], mult*vcat(nans,
                                                   last_val,
                                                   fcast[:forecastpseudo][m_forecast.pseudo_observables[obs], :])[end-84:end], label = "R1. Psi_pi = 300", color = :blue, linestyle = :dash)



            m_forecast = Model1002("ss13")
            m_forecast <= Setting(:forecast_horizons, 32)

            # REGIME 1 PARAMETER VALUE zeta_p with eps exercise (Purple)
            if e == "_eps"
                fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data, param_key = :ζ_p, param_value = mode_use[ind],
                                           param_key2 = :ϵ_p, param_value2 = eps_val)
                nans = repeat([NaN], length(fcast[:histpseudo][m_forecast.pseudo_observables[obs], :]))
                plot!(plots[:plt], dates[end-84:end], mult*vcat(nans,
                                                   fcast[:forecastpseudo][m_forecast.pseudo_observables[obs], :])[end-84:end], label = "R1, eps = $(round(eps_val, digits = 2))", color = :purple)
                plot!(plots[:ait], dates[end-84:end], mult*vcat(nans,
                                                   fcast[:forecastpseudo][m_forecast.pseudo_observables[obs], :])[end-84:end], label = "R1, eps = $(round(eps_val, digits = 2))", color = :purple)
            end

            # PLT. R1 value
            m_forecast = Model1002("ss13")
            m_forecast <= Setting(:forecast_horizons, 32)
            m_forecast <= Setting(:alternative_policy, AltPolicy(:plt, plt_eqcond, plt_solve, forecast_init = plt_forecast_init))
            m_forecast <= Setting(:pgap_value, gap_value)
            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data, param_key = :ζ_p, param_value = mode_use[ind])
            nans = repeat([NaN], length(fcast[:histpseudo][m_forecast.pseudo_observables[obs], :]))
            plot!(plots[:plt], dates[end-84:end], mult*vcat(nans,
                                                   last_val,
                                                   fcast[:forecastpseudo][m_forecast.pseudo_observables[obs], :])[end-84:end], label = "Pre 1990 Slope, PLT", color = :blue, linestyle = :dot)

            # PLT. Eps
            if e == "_eps"
                m_forecast = Model1002("ss13")
                m_forecast <= Setting(:forecast_horizons, 32)
                m_forecast <= Setting(:alternative_policy, AltPolicy(:plt, plt_eqcond, plt_solve, forecast_init = plt_forecast_init))
                m_forecast <= Setting(:pgap_value, gap_value)
                fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data, param_key = :ζ_p, param_value = mode_use[ind],
                                           param_key2 = :ϵ_p, param_value2 = eps_val)
                nans = repeat([NaN], length(fcast[:histpseudo][m_forecast.pseudo_observables[obs], :]))
                plot!(plots[:plt], dates[end-84:end], mult*vcat(nans,
                                                   last_val,
                                                   fcast[:forecastpseudo][m_forecast.pseudo_observables[obs], :])[end-84:end], label = "PLT eps = $(round(eps_val, digits = 2))", color = :purple, linestyle = :dot)
            end

            # AIT. R1 Value
            m_forecast = Model1002("ss13")
            m_forecast <= Setting(:forecast_horizons, 32)
            m_forecast <= Setting(:alternative_policy, AltPolicy(:ait, ait_eqcond, ait_solve, forecast_init = ait_forecast_init))
            m_forecast <= Setting(:pgap_value, gap_value)
            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data, param_key = :ζ_p, param_value = mode_use[ind])
            plot!(plots[:ait], dates[end-84:end], mult*vcat(nans,
                                                   last_val,
                                                   fcast[:forecastpseudo][m_forecast.pseudo_observables[obs], :])[end-84:end], label = "Pre 1990 Slope, AIT", color = :blue, linestyle = :dash)

            # AIT. Eps
            if e == "_eps"
                m_forecast = Model1002("ss13")
                m_forecast <= Setting(:forecast_horizons, 32)
                m_forecast <= Setting(:alternative_policy, AltPolicy(:ait, ait_eqcond, ait_solve, forecast_init = ait_forecast_init))
                m_forecast <= Setting(:pgap_value, gap_value)
                fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data, param_key = :ζ_p, param_value = mode_use[ind],
                                           param_key2 = :ϵ_p, param_value2 = eps_val)
                plot!(plots[:ait], dates[end-84:end], mult*vcat(nans,
                                                   last_val,
                                                   fcast[:forecastpseudo][m_forecast.pseudo_observables[obs], :])[end-84:end], label = "AIT eps = $(round(eps_val, digits = 2))", color = :purple, linestyle = :dash)
            end


            # PLT. R2 Value
            m_forecast  = Model1002("ss13")
            m_forecast <= Setting(:forecast_horizons, 32)
            m_forecast <= Setting(:alternative_policy, AltPolicy(:plt, plt_eqcond, plt_solve, forecast_init = plt_forecast_init))
            m_forecast <= Setting(:pgap_value, gap_value)
            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data)
            p_plt = plot!(plots[:plt], dates[end-84:end], mult*vcat(nans,
                                                           last_val,
                                                           fcast[:forecastpseudo][m_forecast.pseudo_observables[obs], :])[end-84:end],
                          label = "Post 1990 Slope, PLT", color = :red, linestyle = :dot)

            # AIT. R2 Value
            m_forecast = Model1002("ss13")
            m_forecast <= Setting(:forecast_horizons, 32)
            m_forecast <= Setting(:alternative_policy, AltPolicy(:ait, ait_eqcond, ait_solve, forecast_init = ait_forecast_init))
            m_forecast <= Setting(:pgap_value, gap_value)
            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data)
            p_ait = plot!(plots[:ait], dates[end-84:end], mult*vcat(nans,
                                                           last_val,
                                                           fcast[:forecastpseudo][m_forecast.pseudo_observables[obs], :])[end-84:end],
                          label = "Post 1990 Slope, AIT", color = :red, linestyle = :dash)

            savefig(plots[:ait], "$(savepath)/output_data/m1002/ss10/forecast/figures/$(DSGE.detexify(obs))_forecast_$(yr)_$(mode)_$(gap_value)_ait$(e).pdf")
        else
            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars, mode_use[1:end .!=ind2], data)
            sys = compute_system(m_forecast)
            irf_s, irf_o, irf_p = impulse_responses(m_forecast, sys)
            irf_plot_b_sh = plot(irf_o[m_forecast.observables[obs], :, m_forecast.exogenous_shocks[:b_sh]], title = "$(string(DSGE.detexify(obs))) b shock", label = "R1", left_margin = 30px, color = :blue)

            irf_plot_rm_sh = plot(irf_o[m_forecast.observables[obs], :, m_forecast.exogenous_shocks[:rm_sh]], title = "$(string(DSGE.detexify(obs))) rm shock", label = "R1", left_margin = 30px, color = :blue)

            dates = vcat(plot_start_date:Dates.Month(3):Date(2030,3,31))
            dates = map(x -> "$(year(x))-Q$(quarterofyear(x))", dates)
            dat = data[data[:date] .>= plot_start_date, :]
            for p in keys(plots)
                plots[p] = plot(dates, 4*vcat(dat[:, obs], repeat([NaN], 60)), label = "", color = :black, legend = :bottomright, left_margin = 30px)
                plot!(plots[p], ylims = (-0.1, 3.0))
                plot!(plots[p], dates, 4*vcat(repeat([NaN], size(dat, 1)-1),
                       data[end, obs],
                       fcast[:forecastobs][m_forecast.observables[obs],:]), label = "Pre 1990 Slope",
                      color = :blue, xticks = (1:20:length(dates), dates[1:20:length(dates)]))
            end

            # R1, raise psi pi
            m_forecast = Model1002("ss13")
            m_forecast <= Setting(:forecast_horizons, 32)
            m_forecast <= Setting(:alternative_policy, AltPolicy(:zero_pi_star, eqcond_zero_pi_star,
                                                                 solve_zero_pi_star,
                                                                 forecast_init = forecast_init_zero_pi_star))
            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           mode_use[1:end .!=ind2], data,
                                           param_key = :ψ1, param_value = 300.)
            plot!(plots[:psi], dates, 4*vcat(repeat([NaN], size(dat, 1)-1),
                       data[end, obs],
                       fcast[:forecastobs][m_forecast.observables[obs],:]), label = "R1, psi_pi = 300", color = :blue,
                  linestyle = :dash)
            m_forecast = Model1002("ss13")
            m_forecast <= Setting(:forecast_horizons, 32)


            # R1 zeta_p with eps exercise (purple)
            if e == "_eps"
                m = Model1002("ss13")
                m_forecast <= Setting(:forecast_horizons, 32)

                fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data, param_key = :ζ_p, param_value = mode_use[ind],
                                           param_key2 = :ϵ_p, param_value2 = eps_val)
                 @show fcast[:forecastobs][m_forecast.observables[obs], :]
                 plot!(plots[:plt], dates, 4*vcat(repeat([NaN], size(dat, 1)-1),
                                                   dat[end, obs],
                                                   fcast[:forecastobs][m_forecast.observables[obs],:]),
                  label = "R1, eps = $(round(eps_val, digits = 2))", color = :purple)
                 plot!(plots[:ait], dates, 4*vcat(repeat([NaN], size(dat, 1)-1),
                                                   dat[end, obs],
                                                   fcast[:forecastobs][m_forecast.observables[obs],:]),
                  label = "R1, eps = $(round(eps_val, digits = 2))", color = :purple)
            end

            # PLT. R1 Value
            m_forecast <= Setting(:alternative_policy, AltPolicy(:plt, plt_eqcond, plt_solve, forecast_init = plt_forecast_init))
            m_forecast <= Setting(:forecast_horizons, 32)
            m_forecast <= Setting(:pgap_value, gap_value)

            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars, mode_use[1:end .!=ind2], data)
            plot!(plots[:plt], dates, 4*vcat(repeat([NaN], size(dat, 1)-1),
                                 dat[end, obs],
                                 fcast[:forecastobs][m_forecast.observables[obs],:]), label = "Pre 1990 Slope, PLT",
                                 color = :blue, linestyle = :dot)

            # PLT. Eps
            if e == "_eps"
                 m_forecast = Model1002("ss13")
                 m_forecast <= Setting(:forecast_horizons, 32)
                 m_forecast <= Setting(:alternative_policy, AltPolicy(:plt, plt_eqcond, plt_solve, forecast_init = plt_forecast_init))
                 m_forecast <= Setting(:pgap_value, gap_value)

                 fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars, mode_use[1:end .!=ind2], data,
                                           param_key = :ϵ_p, param_value = eps_val)

                 plot!(plots[:plt], dates, 4*vcat(repeat([NaN], size(dat, 1)-1),
                                 dat[end, obs],
                                 fcast[:forecastobs][m_forecast.observables[obs],:]), label = "PLT Eps",
                                 color = :purple, linestyle = :dot)
            end

            # AIT. R1 Value
            m_forecast = Model1002("ss13")
            m_forecast <= Setting(:forecast_horizons, 32)
            m_forecast <= Setting(:alternative_policy, AltPolicy(:ait, ait_eqcond, ait_solve, forecast_init = ait_forecast_init))
            m_forecast <= Setting(:pgap_value, gap_value)

            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars, mode_use[1:end .!=ind2], data)
            plot!(plots[:ait], dates, 4*vcat(repeat([NaN], size(dat, 1)-1),
                                 dat[end, obs],
                                 fcast[:forecastobs][m_forecast.observables[obs],:]), label = "Pre 1990 Slope, AIT",
                                 color = :blue, linestyle = :dash)
            # AIT. Eps
            if e == "_eps"
                 m_forecast = Model1002("ss13")
                 m_forecast <= Setting(:forecast_horizons, 32)
                 m_forecast <= Setting(:alternative_policy, AltPolicy(:ait, ait_eqcond, ait_solve, forecast_init = ait_forecast_init))
                 m_forecast <= Setting(:pgap_value, gap_value)

                 fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars, mode_use[1:end .!=ind2], data,
                                           param_key = :ϵ_p, param_value = eps_val)
                 plot!(plots[:ait], dates, 4*vcat(repeat([NaN], size(dat, 1)-1),
                                 dat[end, obs],
                                 fcast[:forecastobs][m_forecast.observables[obs],:]), label = "AIT Eps",
                                 color = :purple, linestyle = :dash)
            end

            m_forecast = Model1002("ss13")
            m_forecast <= Setting(:forecast_horizons, 32)
            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data)

            for p in keys(plots)
                plot!(plots[p], dates, 4*vcat(repeat([NaN], size(dat, 1)-1),
                              dat[end, obs],
                              fcast[:forecastobs][m_forecast.observables[obs],:]), label = "Post 1990 Slope", color = :red)
                @show fcast[:forecastobs][m_forecast.observables[obs], :]

            end


            sys = compute_system(m_forecast)
            irf_s, irf_o, irf_p = impulse_responses(m_forecast, sys)
            plot!(irf_plot_b_sh, irf_o[m_forecast.observables[obs], :, m_forecast.exogenous_shocks[:b_sh]], title = "$(string(DSGE.detexify(obs))) b shock", label = "R2", color = :red)

            plot!(irf_plot_rm_sh, irf_o[m_forecast.observables[obs], :, m_forecast.exogenous_shocks[:rm_sh]], title = "$(string(DSGE.detexify(obs))) rm shock", label = "R2", color = :red)


            # R2, raise psi pi
            m_forecast = Model1002("ss13")
            m_forecast <= Setting(:forecast_horizons, 32)
            m_forecast <= Setting(:alternative_policy, AltPolicy(:zero_pi_star, eqcond_zero_pi_star,
                                                                 solve_zero_pi_star,
                                                                 forecast_init = forecast_init_zero_pi_star))
            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data, param_key = :ψ1, param_value = 300.)

            plot!(plots[:psi], dates, 4*vcat(repeat([NaN], size(dat, 1)-1),
                              dat[end, obs],
                              fcast[:forecastobs][m_forecast.observables[obs],:]), label = "R2, psi_pi = 300",
                  color = :red,
                  linestyle = :dash)


            # PLT. R2 Value
            m_forecast = Model1002("ss13")
            m_forecast <= Setting(:forecast_horizons, 32)
            m_forecast <= Setting(:alternative_policy, AltPolicy(:plt, plt_eqcond, plt_solve, forecast_init = plt_forecast_init))
            m_forecast <= Setting(:pgap_value, gap_value)

            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data)
            plot!(plots[:plt], dates, 4*vcat(repeat([NaN], size(dat, 1)-1),
                                   dat[end, obs],
                                   fcast[:forecastobs][m_forecast.observables[obs],:]), label = "Post 1990 Slope, PLT", color = :red, linestyle = :dot)

            # AIT. R2 Value
            m_forecast = Model1002("ss13")
            m_forecast <= Setting(:forecast_horizons, 32)
            m_forecast <= Setting(:alternative_policy, AltPolicy(:ait, ait_eqcond, ait_solve, forecast_init = ait_forecast_init))
            m_forecast <= Setting(:pgap_value, gap_value)

            fcast = DSGE.forecast_one_draw(m_forecast, :mode, :none, output_vars,
                                           vcat(mode_use[1:(ind-1)], mode_use[ind2],
                                                mode_use[(ind+1):(ind2-1)], mode_use[(ind2+1):end]),
                                           data)
            plot!(plots[:ait], dates, 4*vcat(repeat([NaN], size(dat, 1)-1),
                                   dat[end, obs],
                                   fcast[:forecastobs][m_forecast.observables[obs],:]), label = "Post 1990 Slope, AIT", color = :red, linestyle = :dash)

            savefig(plots[:ait], "$(savepath)/output_data/m1002/ss10/forecast/figures/$(DSGE.detexify(obs))_forecast_$(yr)_$(mode)_$(gap_value)_ait$(e).pdf")
        end
        end
    end
end
