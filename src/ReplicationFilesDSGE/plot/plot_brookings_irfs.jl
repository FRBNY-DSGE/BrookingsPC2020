function plot_brookings_irfs(m1::Union{AbstractDSGEModel,AbstractVARModel},
                             m2::Union{AbstractDSGEModel,AbstractVARModel},
                             vars::Vector{Symbol}, class::Symbol,
                             irf_type::Symbol,
                             input_type::Symbol, cond_type::Symbol;
                             irf_product::Symbol = :irf,
                             observables::Vector{Symbol} = Vector{Symbol}(undef,0),
                             forecast_string1::String = "",
                             forecast_string2::String = "",
                             titles::Vector{String} = Vector{String}(undef,0),
                             legend_names::Vector{String} = ["", ""],
                             color1::Union{Symbol,AbstractRGB} = RGB(55. / 255., 126. / 255., 184. / 255.), # :blue
                             color2::Union{Symbol,AbstractRGB} = RGB(.8941, .1020, .1098), # :red
                             alpha1::Float64 = 0.5,
                             alpha2::Float64 = 0.5,
                             density_bands::Vector{String} = ["90.0%"],
                             plotroot::String = figurespath(m1, "forecast"),
                             fileformat::String =
                             string(DSGE.plot_extension()),
                             ylims_dict::Dict{Symbol,Tuple{Float64,Float64}} =
                             Dict{Symbol,Tuple{Float64,Float64}}(),
                             addl_text::String = "",
                             shocks::Vector{Symbol} = Vector{Symbol}(undef,0),
                             lags::Int = 4, obs_shock::Symbol = :none, flip_shocks::Bool = false,
                             which_model::Int = 1, compute_mode::Bool = false,
                             ytickfont::Int = 12, xtickfont::Int = 12,
                             ylabels_dict::Dict{Symbol, String} = Dict{Symbol, String}(),
                             xlabels_dict::Dict{Symbol, String} = Dict{Symbol, String}(),
                             verbose::Symbol = :low) where {S<:Real}
    if !(irf_type in [:cholesky, :maxBC, :max_business_cycle_variance, :cholesky_long_run,
                      :choleskyLR])
        error("Input irf_type cannot be $irf_type. It must be one of " *
              ":cholesky, :maxBC, :max_buiness_cycle_variance, :cholesky_long_run,"
              * " or :choleskyLR")
    end

    # Read in MeansBands
    if (irf_product == :dsgevarirf || irf_product == :dsgevarrotationirf ||
        irf_product == :dsgevarlambdairf)
        var_names = isa(which_model == 1 ? m1 : m2, AbstractDSGEModel) ?
            "_" * join(string.(map(x -> DSGE.detexify(x), observables)), "_") : ""
        mb_name = Symbol(irf_product, :obs, Symbol(var_names), :_, irf_type, addl_text)
    else
        mb_name = Symbol(irf_product, class, :_, irf_type)
    end
    mb1 = read_mb(m1, input_type, cond_type, mb_name;
                  forecast_string = forecast_string1)
    mb2 = read_mb(m2, input_type, cond_type, mb_name;
                  forecast_string = forecast_string2)

    # Compute modal IRFs
    if !isempty(observables) && !isempty(shocks) && compute_mode
        DSGE.update!(m1, load_draws(m1, :mode))
        DSGE.update!(m2, load_draws(m2, :mode))
        mode_irf1 = impulse_responses(m1, reshape(map(x -> x.value, DSGE.get_parameters(m1)),
                                                  1, n_parameters(m1)),
                                      :mode, irf_type,
                                      lags, observables, shocks,
                                      findfirst(observables .== obs_shock);
                                      parallel = false, flip_shocks = flip_shocks)
        mode_irf2 = impulse_responses(m2, reshape(map(x -> x.value, DSGE.get_parameters(m2)),
                                                  1, n_parameters(m2)),
                                      :mode, irf_type,
                                      lags, observables, shocks,
                                      findfirst(observables .== obs_shock);
                                      parallel = false, flip_shocks = flip_shocks)
    end

    # Get titles
    if isempty(titles)
        detexify_title = typeof(Plots.backend()) == Plots.GRBackend
        title_vars = Vector{Symbol}(undef,0)
        for var in vars
            if var in observables
                push!(title_vars, var)
            end
        end
        if which_model == 1
            titles = map(var ->
                         DSGE.describe_series(m1, var, class, detexify = detexify_title),
                         title_vars)
        elseif which_model == 2
            titles = map(var ->
                         DSGE.describe_series(m2, var, class, detexify = detexify_title),
                         title_vars)
        end
        if length(title_vars) < length(vars)
            push!(titles, [string(DSGE.detexify(var)) for var in setdiff(vars, title_vars)]...)
        end
    end

    # Loop and plot
    outplots = OrderedDict{Symbol,Plots.Plot}()
    for (var, title) in zip(vars, titles)
        outplots[var] = Plots.plot()
        h = size(mb1.means,1)
        mb1.means[!,:date] = collect(1:h)
        mb2.means[!,:date] = collect(1:h)
        if isa(color1, AbstractRGB)
            band_color1 = RGBA(color1.r, color1.g, color1.b, alpha1)
        end
        if isa(color2, AbstractRGB)
            band_color2 = RGBA(color2.r, color2.g, color2.b, alpha2)
        end
        tmp = brookingsirf!(var, mb1; mean_color = color1, #title = title,
                           bands_color = band_color1, band_alpha = alpha1,
                           meanlabel = legend_names[1],
                           density_bands = density_bands)
        brookingsirf!(var, mb2; mean_color = color2, bands_color = band_color2,
                     band_alpha = alpha2, meanlabel = legend_names[2],
                     density_bands = density_bands)
        if !isempty(observables) && !isempty(shocks) && var in observables && compute_mode
            mode_i = findfirst(observables .== var)
            plot!(1:h, mode_irf1[mode_i,:], color = color1,
                  linestyle = :dash, label = legend_names[1] * " (Mode)")
            plot!(1:h, mode_irf2[mode_i,:], color = color2,
                  linestyle = :dash, label = legend_names[2] * " (Mode)")
            if xtickfont > 0
                plot!(tmp, xtickfont = font(12))
            end
            if ytickfont > 0
                plot!(tmp, ytickfont = font(12))
            end
        end

        outplots[var] = tmp
    end

    # Change ylims of named plots
    for (k,v) in ylims_dict # if empty, this loop never happens
        if haskey(outplots, k)
            plot(outplots[k])
            ylims!(v)
        end
    end

    # Change ylabel of named plots
    for (k,v) in ylabels_dict # if empty, this loop never happens
        if haskey(outplots, k)
            outplots[k] = plot(outplots[k])
            ylabel!(v)
        end
    end

    # Change xlabel of named plots
    for (k,v) in xlabels_dict # if empty, this loop never happens
        if haskey(outplots, k)
            outplots[k] = plot(outplots[k])
            xlabel!(v)
        end
    end

    # Save plots
    if !isempty(plotroot)
        if irf_product == :dsgevarirf
            tail = isa(which_model == 1 ? m1 : m2, AbstractDSGEModel) ?
                Symbol(irf_type, "_variables", var_names) :
                Symbol(irf_type)
        else
            tail = irf_type
        end
        for var in vars
            try
                if which_model == 1
                    filename = get_forecast_filename(plotroot, filestring_base(m1),
                                                     input_type, cond_type,
                                                     Symbol("brookings", addl_text, "_",  string(irf_product), "_",
                                                            DSGE.detexify(var), "_",
                                                            tail),
                                                     forecast_string = forecast_string1 *
                                                     forecast_string2,
                                                     fileformat = fileformat)
                elseif which_model == 2
                    filename = get_forecast_filename(plotroot, filestring_base(m2),
                                                     input_type, cond_type,
                                                     Symbol("brookings", addl_text, "_",  string(irf_product), "_",
                                                            DSGE.detexify(var), "_",
                                                            tail),
                                                     forecast_string = forecast_string1 *
                                                     forecast_string2,
                                                     fileformat = fileformat)
                end
                DSGE.save_plot(outplots[var], filename, verbose = verbose)
            catch err
                if isa(err, SystemError) && irf_product == :dsgevarirf &&
                    isa(which_model == 1 ? m1 : m2, AbstractDSGEModel)
                    newvar_names = join(string.(map(x -> DSGE.detexify(x), observables)), "_")
                    @show newvar_names
                    newvar_names = Symbol(replace(newvar_names, "obs_" => ""))
                    newtail = Symbol(irf_type, "_variables_", newvar_names)
                    if which_model == 1
                        filename = get_forecast_filename(plotroot, filestring_base(m1),
                                                         input_type, cond_type,
                                                         Symbol("brookings", addl_text, "_",  string(irf_product), "_",
                                                                DSGE.detexify(var), "_",
                                                                newtail),
                                                         forecast_string = "",
                                                         fileformat = fileformat)
                    else
                        filename = get_forecast_filename(plotroot, filestring_base(m2),
                                                         input_type, cond_type,
                                                         Symbol("brookings", addl_text, "_",  string(irf_product), "_",
                                                                DSGE.detexify(var), "_",
                                                                newtail),
                                                         forecast_string = "",
                                                         fileformat = fileformat)
                    end
                    DSGE.save_plot(outplots[var], filename, verbose = verbose)
                else
                    throw(err)
                end
            end
        end
    end

    return outplots
end

function plot_brookings_irfs(m::AbstractDSGEModel, vars::Vector{Symbol}, class::Symbol,
                             irf_type::Symbol, input_type::Symbol, cond_type::Symbol;
                             forecast_string::String = "",
                             titles::Vector{String} = Vector{String}(undef,0),
                             legend_name::String = "",
                             color::Symbol = :blue,
                             alpha::Float64 = 0.5,
                             density_bands::Vector{String} = ["90.0%"],
                             plotroot::String = figurespath(m, "forecast"),
                             fileformat::String =
                             string(DSGE.plot_extension()),
                             ylims_dict::Dict{Symbol,Tuple{Float64,Float64}} =
                             Dict{Symbol,Tuple{Float64,Float64}}(),
                             verbose::Symbol = :low) where {S<:Real}
    if !(irf_type in [:cholesky, :maxBC, :max_business_cycle_variance])
        error("Input irf_type cannot be $irf_type. It must be one of " *
              ":cholesky, :maxBC, or :max_buiness_cycle_variance")
    end


    # Read in MeansBands
    mb = read_mb(m, input_type, cond_type, Symbol(:irf, class, :_, irf_type),
                  forecast_string = forecast_string)

    # Get titles
    if isempty(titles)
        detexify_title = typeof(Plots.backend()) == Plots.GRBackend
        titles = map(var ->
                     DSGE.describe_series(m, var, class, detexify = detexify_title), vars)
    end

    # Loop and plot
    outplots = OrderedDict{Symbol,Plots.Plot}()
    for (var, title) in zip(vars, titles)
        outplots[var] = Plots.plot()
        h = size(mb.means,1)
        mb.means[!,:date] = collect(1:h)
        tmp = brookingsirf!(var, mb; mean_color = color, title = title,
                           bands_color = color, band_alpha = alpha,
                           meanlabel = legend_name,
                           density_bands = density_bands)
        outplots[var] = tmp
    end

    # Change ylims of named plots
    for (k,v) in ylims_dict # if empty, this loop never happens
        if haskey(outplots, k)
            plot(outplots[k])
            ylims!(v)
        end
    end

    # Save plots
    if !isempty(plotroot)
        var_names
        for var in vars
            filename = get_forecast_filename(plotroot, filestring_base(m),
                                             input_type, cond_type,
                                             Symbol("brookings_irf",
                                                    "_",
                                                    DSGE.detexify(var), "_",
                                                    irf_type, "_",
                                                    "varobs_", ),
                                             forecast_string = forecast_string,
                                             fileformat = fileformat)
            DSGE.save_plot(outplots[var], filename, verbose = verbose)
        end
    end

    return outplots
end

@userplot BrookingsIrf

brookingsirf

@recipe function f(brookingsirf::BrookingsIrf;
                   mean_color = :blue,
                   bands_color = :blue,
                   density_bands = which_density_bands(brookingsirf.args[2], uniquify = true),
                   band_alpha = 0.1,
                   meanlabel = "")
    gr(display_type = :inline)
    var, mb = brookingsirf.args
    quarters_ahead = collect(1:size(mb.means,1))

    # Bands
    for pct in density_bands
        @series begin
            fillcolor := bands_color
            fillalpha := band_alpha
            linealpha := 0
            label     := meanlabel #""
            left_margin := 40px
            bottom_margin := 40px
            fillrange := mb.bands[var][!,Symbol(pct, " UB")]
            quarters_ahead, mb.bands[var][!,Symbol(pct, " LB")]
        end
    end

    # Mean
    @series begin
        linewidth := 2
        linecolor := mean_color
        label     := "" #meanlabel
        alpha     := 1.0
        left_margin := 40px
        bottom_margin := 40px
        quarters_ahead, mb.means[!, var]
    end
end

function plot_brookings_irfs(m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                             obs1::AbstractArray{S}, pseudo1::AbstractArray{S},
                             obs2::AbstractArray{S}, pseudo2::AbstractArray{S},
                             input_type::Symbol, cond_type::Symbol,
                             output_vars::Vector{Symbol};
                             observables_order::Vector{<:Int} =
                             collect(values(m1.observables)),
                             append_stationary_state::Bool = true,
                             forecast_string1::String = "",
                             forecast_string2::String = "",
                             color1::Symbol = :blue,
                             color2::Symbol = :red,
                             legend_names::Vector{String} = ["", ""],
                             plotroot::String = figurespath(m1, "forecast"),
                             fileformat::String =
                             string(DSGE.plot_extension()),
                             verbose::Symbol = :low) where {S<:Real}

    outplots_obs    = Dict{Symbol,Plots.Plot}(var => plot() for var in keys(m1.observables))
    outplots_pseudo = Dict{Symbol,Plots.Plot}(var => plot()
                                              for var in keys(m1.pseudo_observables))

    if append_stationary_state
        obs1 = hcat(zeros(size(obs1,1)), obs1)
        pseudo1 = hcat(zeros(size(pseudo1,1)), pseudo1)
        obs2 = hcat(zeros(size(obs2,1)), obs2)
        pseudo2 = hcat(zeros(size(pseudo2,1)), pseudo2)

        # How much to subtract from x values for plotting
        start_x = 1
    else
        start_x = 0
    end

    # Set x values for plotting
    xvar_obs1 = collect(1:size(obs1,2)) .- start_x
    xvar_obs2 = collect(1:size(obs2,2)) .- start_x
    xvar_pseudo1 = collect(1:size(pseudo1,2)) .- start_x
    xvar_pseudo2 = collect(1:size(pseudo2,2)) .- start_x



    # Mode or Full?
    if input_type == :mode
        if :irfobs in output_vars
            T = maximum([size(obs1,2), size(obs2,2)])
            for var in collect(keys(m1.observables))[observables_order]
                tmp = plot(outplots_obs[var])
                plot!(xvar_obs1, obs1[m1.observables[var],:],
                      label = legend_names[1], linewidth = 2, color = color1)
                plot!(xvar_obs2, obs2[m2.observables[var],:],
                      label = legend_names[2], linewidth = 2, color = color2,
                      xlim = [0, T])
              #  title!(string(DSGE.detexify(var)))
                xticks!(collect(1:4:T), string.(1:4:(T-1)))
                outplots_obs[var] = tmp
            end
        end

        if :irfpseudo in output_vars
            T = maximum([size(pseudo1,2), size(pseudo2,2)])
            for var in keys(m1.pseudo_observables)
                tmp = plot(outplots_pseudo[var])
                plot!(xvar_pseudo1, pseudo1[m1.pseudo_observables[var],:],
                      label = legend_names[1], linewidth = 2, color = color1)
                plot!(xvar_pseudo2, pseudo2[m2.pseudo_observables[var],:],
                      label = legend_names[2], linewidth = 2, color = color2,
                      xlim = [0, T])
             #   title!(string(DSGE.detexify(var)))
                xticks!(1 .+ collect(1:4:T), string.(1:4:T-1))
                outplots_pseudo[var] = tmp
            end
        end

    else
        if :irfobs in output_vars
            T = maximum([size(obs1,2), size(obs2,2)])
            for var in collect(keys(m1.observables))[observables_order]
                tmp = plot(outplots_obs[var])
                plot!(xvar_obs1, obs1[m1.observables[var],:],
                      label = legend_names[1], linewidth = 2, color = color1)
                plot!(xvar_obs2, obs2[m2.observables[var],:],
                      label = legend_names[2], linewidth = 2, color = color2,
                      xlim = [0, T])
            #    title!(string(DSGE.detexify(var)))
                xticks!(collect(1:4:T), string.(1:4:(T-1)))
                outplots_obs[var] = tmp
            end
        end

        if :irfpseudo in output_vars
            T = maximum([size(pseudo1,2), size(pseudo2,2)])
            for var in keys(m1.pseudo_observables)
                tmp = plot(outplots_pseudo[var])
                plot!(xvar_pseudo1, pseudo1[m1.pseudo_observables[var],:],
                      label = legend_names[1], linewidth = 2, color = color1)
                plot!(xvar_pseudo2, pseudo2[m2.pseudo_observables[var],:],
                      label = legend_names[2], linewidth = 2, color = color2,
                      xlim = [0, T])
             #   title!(string(DSGE.detexify(var)))
                xticks!(1 .+ collect(1:4:T), string.(1:4:T-1))
                outplots_pseudo[var] = tmp
            end
        end


    end
    if :irfstates in output_vars
        @warn ":irfstates has not been implemented yet. Skipping these plots..."
    end

    if !isempty(plotroot)
        for var in keys(m1.observables)
            obs_file = get_forecast_filename(plotroot, filestring_base(m1),
                                             input_type, cond_type,
                                             Symbol("brookings_irf", "_",
                                                    DSGE.detexify(var)),
                                             forecast_string = forecast_string1 *
                                             forecast_string2,
                                             fileformat = fileformat)
            DSGE.save_plot(outplots_obs[var], obs_file, verbose = verbose)
        end
        for var in keys(m1.pseudo_observables)
            pseudo_file = get_forecast_filename(plotroot, filestring_base(m1),
                                             input_type, cond_type,
                                             Symbol("brookings_irf", "_",
                                                    DSGE.detexify(var)),
                                             forecast_string = forecast_string1 *
                                             forecast_string2,
                                             fileformat = fileformat)
            DSGE.save_plot(outplots_pseudo[var], pseudo_file, verbose = verbose)
        end
    end

    return outplots_obs, outplots_pseudo
end
