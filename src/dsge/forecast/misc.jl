# Correct log labor productivity and log wages by accumulating z_t
function correct_integ_series!(m::AbstractDSGEModel,
                                   pseudo::Matrix{T}) where {T<:AbstractFloat}
    lpi = m.pseudo_observables[:laborproductivity]
    wti = m.pseudo_observables[:Wages]
    zti = m.pseudo_observables[:z_t]
    z_t_cum = cumsum(pseudo[zti,:])
    pseudo[lpi,:] += z_t_cum
    pseudo[wti,:] += z_t_cum
    return pseudo
end

# Compute MeansBands from saved irf output
function compute_correct_integ_series_mb(m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                                         input_type1::Symbol, input_type2::Symbol,
                                         cond_type1::Symbol, cond_type2::Symbol,
                                         shocks::Vector{Symbol};
                                         forecast_string1::String = "",
                                         forecast_string2::String = "",
                                         addl_text::String = "",
                                         which_model::Int = 1,
                                         density_bands::Vector{Float64} = [.5, .6, .7, .8, .9])
    integ_serieses = [:Wages, :laborproductivity]
    for integ_series in integ_serieses
        for shock in shocks
            println("Correcting IRF of $integ_series to shock $shock")
            irf1, _ = read_forecast_output(m1, input_type1, cond_type1,
                                           :irfpseudo, integ_series,
                                           Nullable(shock), forecast_string = forecast_string1)
            irf2, _ = read_forecast_output(m2, input_type2, cond_type2,
                                           :irfpseudo, integ_series,
                                           Nullable(shock), forecast_string = forecast_string2)
            zt1, _ = read_forecast_output(m1, input_type1, cond_type1,
                                          :irfpseudo, :z_t,
                                          Nullable(shock), forecast_string = forecast_string1)
            zt2, _ = read_forecast_output(m2, input_type2, cond_type2,
                                          :irfpseudo, :z_t,
                                          Nullable(shock), forecast_string = forecast_string2)

            # Correct labor productivity
            z_t1_cum = cumsum(zt1, dims = 2)
            z_t2_cum = cumsum(zt2, dims = 2)
            irf1 += z_t1_cum
            irf2 += z_t2_cum

            # Compute means and bands
            means1_data = vec(mean(irf1, dims = 1))
            means2_data = vec(mean(irf2, dims = 1))
            bands1_data = DSGE.find_density_bands(irf1, density_bands, minimize = true)
            bands2_data = DSGE.find_density_bands(irf2, density_bands, minimize = true)
            means1 = DataFrame()
            bands1 = Dict{Symbol,DataFrame}()
            means2 = DataFrame()
            bands2 = Dict{Symbol,DataFrame}()
            means1[!, integ_series] = means1_data
            means2[!, integ_series] = means2_data
            bands1[integ_series] = bands1_data
            bands2[integ_series] = bands2_data


            metadata1 = Dict{Symbol,Any}()
            metadata2 = Dict{Symbol,Any}()
            names1 = keys(m1.pseudo_observables)
            names2 = keys(m2.pseudo_observables)
            metadata1[:para] = input_type1
            metadata1[:cond_type] = cond_type1
            metadata2[:para] = input_type1
            metadata2[:cond_type] = cond_type1
            metadata1[:product] = :irf
            metadata2[:product] = :irf
            metadata1[:class] = :pseudo
            metadata2[:class] = :pseudo
            metadata1[:indices] = OrderedDict{Symbol,Int}(integ_series => 1)
            metadata2[:indices] = OrderedDict{Symbol,Int}(integ_series => 1)

            mb1 = MeansBands(metadata1, means1, bands1)
            mb2 = MeansBands(metadata2, means2, bands2)

            fp1 = get_meansbands_output_file(m1, input_type1, cond_type1,
                                             Symbol(:irf, :pseudo, :_, shock, :_,
                                                    integ_series,
                                                    Symbol(addl_text)),
                                             forecast_string = forecast_string1)
            fp2 = get_meansbands_output_file(m2, input_type2, cond_type2,
                                             Symbol(:irf, :pseudo, :_, shock, :_,
                                                    integ_series,
                                                    Symbol(addl_text)),
                                             forecast_string = forecast_string2)

            dirpath1 = dirname(fp1)
            dirpath2 = dirname(fp2)
            isdir(dirpath1) || mkpath(dirpath1)
            isdir(dirpath2) || mkpath(dirpath2)
            JLD2.jldopen(fp1, true, true, true, IOStream) do file
                write(file, "mb", mb1)
            end
            JLD2.jldopen(fp2, true, true, true, IOStream) do file
                write(file, "mb", mb2)
            end
        end
    end
end

# plot wages and labor productivity results
function plot_correct_integ_series_irfs(m1::AbstractDSGEModel, m2::AbstractDSGEModel,
                                        shocks::Vector{Symbol},
                                        input_type::Symbol, cond_type::Symbol;
                                        forecast_string1::String = "",
                                        forecast_string2::String = "",
                                        titles::Vector{String} = Vector{String}(undef,0),
                                        legend_names::Vector{String} = ["", ""],
                                        color1::Symbol = :blue,
                                        color2::Symbol = :red,
                                        alpha1::Float64 = 0.5,
                                        alpha2::Float64 = 0.5,
                                        density_bands::Vector{String} = ["90.0%"],
                                        plotroot::String = figurespath(m1, "forecast"),
                                        fileformat::String =
                                        string(DSGE.plot_extension()),
                                        ylims_dict::Dict{Symbol,Tuple{Float64,Float64}} =
                                        Dict{Symbol,Tuple{Float64,Float64}}(),
                                        addl_text::String = "",
                                        verbose::Symbol = :low) where {S<:Real}
    class = :pseudo
    integ_serieses = [:Wages, :laborproductivity]
    integ_serieses_outplots = OrderedDict{Symbol,OrderedDict{Symbol,Plots.Plot}}()
    for integ_series in integ_serieses
        # Get titles
        if isempty(titles)
            detexify_title = typeof(Plots.backend()) == Plots.GRBackend
            title = DSGE.describe_series(m1, integ_series, class, detexify = detexify_title)
        end
        outplots = OrderedDict{Symbol,Plots.Plot}()

        for shock in shocks
            # Read in MeansBands
            mb1 = read_mb(m1, input_type, cond_type, Symbol(:irf, class, :_, shock,
                                                            :_, integ_series, Symbol(addl_text)),
                          forecast_string = forecast_string1)
            mb2 = read_mb(m2, input_type, cond_type, Symbol(:irf, class, :_, shock,
                                                            :_, integ_series, Symbol(addl_text)),
                          forecast_string = forecast_string2)

            # Loop and plot
            outplots[shock] = Plots.plot()
            h = size(mb1.means,1)
            mb1.means[!,:date] = collect(1:h)
            mb2.means[!,:date] = collect(1:h)
            tmp = brookingsirf!(integ_series, mb1; mean_color = color1, title = title,
                                bands_color = color1, band_alpha = alpha1,
                                meanlabel = legend_names[1],
                                density_bands = density_bands)
            brookingsirf!(integ_series, mb2; mean_color = color2, bands_color = color2,
                          band_alpha = alpha2, meanlabel = legend_names[2],
                          density_bands = density_bands)
            outplots[shock] = tmp
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
            for shock in shocks
                filename = get_forecast_filename(plotroot, filestring_base(m1),
                                                 input_type, cond_type,
                                                 Symbol("irf", "_",
                                                        DSGE.detexify(shock), "_",
                                                        DSGE.detexify(integ_series),
                                                        addl_text),
                                                 forecast_string =
                                                 forecast_string1 * forecast_string2,
                                                 fileformat = fileformat)
                DSGE.save_plot(outplots[shock], filename, verbose = verbose)
            end
        end
        integ_serieses_outplots[integ_series] = outplots
    end
    return integ_serieses_outplots
end
