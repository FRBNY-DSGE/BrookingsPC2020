using DSGE, Plots, SMC, JLD2, FileIO, Statistics, Dates, ModelConstructors
using DSGEModels, Printf, DataFrames, ColorTypes, Nullables
using Plots.PlotMeasures
using KernelDensity

save_orig = true
fig_type = "svg"

outdir = if save_orig
    "save_orig"
else
    "save"
end

include("../df_to_tex.jl")
include("../util.jl")
m = Model1002()
gr()

fp = dirname(@__FILE__)
vint = "191118"

function add_tex(fid::IOStream, savename::String)
    @printf fid "\\includegraphics[width=.5\\textwidth]{%s} \n" savename
end


fid = open("$fp/../../../docs/dsge_parameters/dsge_parameters_separate_estimations.tex", "w")
mdds = DataFrame(model = String[], prior = String[], period = String[], half = String[], mdd = Float64[])


for model in ["m1002"]
   if model == "smets_wouters_orig"
        @printf fid "\\section{SW Orig -- Half 1, Half 2}"
    else
        @printf fid "\\section{%s -- Half 1, Half2}" model
    end
    for pr in ["diffuse", "standard"]
        @printf fid "\\subsection{%s}" pr
    if model == "m1002"
        if pr == "diffuse"
            m = Model1002("ss5")
        elseif pr == "standard"
            m = Model1002("ss10")
        end
    elseif model == "m904"
        if pr == "diffuse"
            m = Model904("ss7")
        elseif pr == "standard"
            m = Model904("ss9")
        end
    elseif model == "m805"
        if pr == "diffuse"
            m = Model805("ss5")
        elseif pr == "standard"
            m = Model805("ss1")
        end
    elseif model == "smets_wouters_orig"
        if pr == "diffuse"
            m = SmetsWoutersOrig("ss3")
        elseif pr == "standard"
            m = SmetsWoutersOrig("ss0")
        end
    end
        standard_spec!(m, vint, fp; fcast_date = Date(2019, 12, 31), dsid = 10021, cdid = 1)
        m <= Setting(:saveroot, "../../../$(outdir)/")
        m <= Setting(:npart, "15000", true, "npart", "number of SMC particles")
        m <= Setting(:use_population_forecast, false)

        inds = findall(x -> x.key in [:ζ_p, :ψ1, :ψ2, :ψ3, :ρ, :ι_w, :ι_p, :ζ_w], m.parameters)

        for period in ["full_incZLB", "full_preZLB"]
            if period == "full_incZLB"
                m <= Setting(:preZLB, false, true, "preZLB", "")
                ZLB1 = "incZLB"
                reg_splits = ["900331"]
                @printf fid "\\subsubsection{Including ZLB}"
            elseif period == "full_preZLB"
                m <= Setting(:preZLB, true, true, "preZLB", "")
                ZLB1 = "preZLB"
                reg_splits = ["840331", "900331"]
                @printf fid "\\subsubsection{Pre ZLB}"
            end
            for reg_split in reg_splits
                @printf fid "\\textbf{Split in %s} \\\\" reg_split
                if reg_split == "840331"
                    m <= Setting(:reg2start, reg_split, true,"reg2start", "")
                elseif reg_split == "900331"
                    m <= Setting(:reg2start, reg_split, true, "reg2start", "")
                end
                ZLB = "$(ZLB1)_$(reg_split)"
                try
                    m <= Setting(:period, "r1", true, "period", "")
                    mdd = marginal_data_density(m)
                    push!(mdds, [model, pr, period, "r1", mdd])
                    draws_r1 = load_draws(m, :full)
                    input_file_name = get_forecast_input_file(m, :full)
                    cloud = load(replace(replace(input_file_name, ".h5" => ".jld2"), "smcsave" => "smc_cloud"), "cloud")
                    postmode_r1 = load_draws(m, :mode)

                    m <= Setting(:period, "r2", true, "period", "")

                    mdd = marginal_data_density(m)
                    push!(mdds, [model, pr, period, "r2", mdd])
                    draws_r2 = load_draws(m, :full)
                    input_file_name = get_forecast_input_file(m, :full)
                    cloud = load(replace(replace(input_file_name, ".h5" => ".jld2"), "smcsave" => "smc_cloud"), "cloud")
                    postmode_r2 = load_draws(m, :mode)
                    figure_save_path = "$fp/../../../$(outdir)/output_data/$(spec(m))/$(subspec(m))/estimate/figures/"
                    @show figure_save_path
                    @show ZLB
                    for ind in inds
                        if  pr == "standard" && period == "full_incZLB" && reg_split == "900331"
                            @show m.parameters[ind].key, postmode_r1[ind], postmode_r2[ind]
                            @show "kappa_p", get_κp(m, postmode_r1), get_κp(m, postmode_r2)
                            @show "kappa_w", get_κw(m, postmode_r1), get_κw(m, postmode_r2)
                        end
                        max1 = maximum(draws_r1[:, ind])
                        max2 = maximum(draws_r2[:, ind])
                        max12 = max(max1, max2)
                        min1 = minimum(draws_r1[:, ind])
                        min2 = minimum(draws_r2[:, ind])
                        min12 = min(min1, min2)
                        padding = (max12 - min12) / 100.
                        start = min12 - padding
                        stop  = max12 + padding

                        if m.parameters[ind].key == :ζ_p
                            start -= .4
                        elseif m.parameters[ind].key == :ζ_w
                            start -= .3
                        elseif m.parameters[ind].key == :ψ3
                            start -= .075
                        elseif m.parameters[ind].key == :ρ
                            stop += .05
                        end
                        param_label = replace(replace(replace(DSGE.detexify(m.parameters[ind].tex_label),
                                                              "\\" => "", ), "{" => ""), "}" => "")
                        p = histogram(1:size(draws_r1, 1), draws_r1[:, ind], color = RGB(55. / 255., 126. / 255., 184. / 255.), alpha = 0.5, label = "Pre 1990",
                                      bins = 100,
                                      normalize = :pdf, xlims = (start, stop),
                                      left_margin=20px)
                        histogram!(p, 1:size(draws_r2, 1), draws_r2[:, ind], color = RGB(.8941, .1020, .1098),
                                   alpha = 0.5, label = "Post 1990", bins = 100, normalize = :pdf,
                                   xlims = (start, stop), left_margin=20px, xtickfont = font(12), ytickfont = font(12))
                        plot!(p, range(start, length = 1000, stop = stop),
                              pdf.(get(m.parameters[ind].prior),
                                   range(start, length = 1000, stop = stop)),
                              label = "Prior", color = :black, normalize = :pdf,
                              xlims = (start, stop))

                        savefig(p, "$(figure_save_path)/$(param_label)_$(ZLB).$(fig_type)")
                        add_tex(fid, "$(fp)/../../../$(outdir)/output_data/$(model)/$(subspec(m))/estimate/figures/$(param_label)_$(ZLB).(fig_type)")
                    end
                    if model == "m1002"
                        kdep, kdew = draw_κ(m)
                        for i in [:κ_p, :κ_w]
                            param_label = string(DSGE.detexify(i))
                            r1draws, r2draws, r1postmode, r2postmode = if i == :κ_p
                                get_κp(m, draws_r1), get_κp(m, draws_r2),
                                get_κp(m, reshape(postmode_r1, 1, DSGE.n_parameters(m))),
                                get_κp(m, reshape(postmode_r2, 1, DSGE.n_parameters(m)))
                            else
                                get_κw(m, draws_r1), get_κw(m, draws_r2),
                                get_κw(m, reshape(postmode_r1, 1, DSGE.n_parameters(m))),
                                get_κw(m, reshape(postmode_r2, 1, DSGE.n_parameters(m)))
                            end
                            p = histogram(1:size(r1draws, 1), r1draws, color = RGB(55. / 255., 126. / 255., 184. / 255.),
                                          alpha = 0.5, label = "Pre 1990",
                                          bins = 100,
                                          normalize = :pdf,
                                          left_margin=20px,
                                          xtickfont = font(12), ytickfont = font(12),
                                          right_margin = 10px)
                            histogram!(p, 1:size(r2draws, 1), r2draws, color = color = RGB(.8941, .1020, .1098), alpha = 0.5,
                                       label = "Post 1990", bins = 100, normalize = :pdf,
                                       left_margin=20px, xtickfont = font(12), ytickfont = font(12),
                                       right_margin = 10px)
                            max1 = maximum(r1draws)
                            max2 = maximum(r2draws)
                            max12 = max(max1, max2)
                            min1 = minimum(r1draws)
                            min2 = minimum(r2draws)
                            min12 = min(min1, min2)
                            padding = (max12 - min12) / 50.
                            start = min12 - padding
                            stop  = max12 + padding
                            if i == :κ_p
                                stop = .1
                            elseif i == :κ_w
                                stop =  .1
                            end

                            if i == :κ_p
                                plot!(xlims = (start, stop))
                               plot!(p, kdep.x, kdep.density, label = "Prior", color = :black,
                                      xlims = (start, stop), normalize = :pdf)
                            else
                                plot!(xlims = (start, stop))
                                plot!(p, kdew.x, kdew.density, label = "Prior", color = :black,
                                      xlims = (start, stop), normalize = :pdf)
                            end
                            savefig(p, "$(figure_save_path)/$(param_label)_$(ZLB).$(fig_type)")
                            add_tex(fid, "$(fp)/../../../$(outdir)/output_data/$(model)/$(subspec(m))/estimate/figures/$(param_label)_$(ZLB).$(fig_type)")
                        end
                    end
                catch err
                    if isa(err, SystemError)
                        nothing
                    else
                        rethrow(err)
                    end
                end
            end
        end
    end
end

close(fid)
