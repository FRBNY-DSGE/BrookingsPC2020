using DSGE, Plots, SMC, JLD2, FileIO, Statistics, Dates, ModelConstructors
using Plots.PlotMeasures, ColorTypes
using DSGEModels, Printf, DataFrames

replicate_exact = true

outer = nothing
include("../df_to_tex.jl")
include("../util.jl")

gr()
GR.inline("pdf")
GR.inline("png")
fp = dirname(@__FILE__)
vint = "191118"

function plot_regime_hist_comp(m::AbstractDSGEModel, draws::Matrix{Float64}, mode::Vector{Float64}, savename::String, ind1::Int, ind2::Int, title::String)
    @show mean(draws[:, ind1])
    @show mean(draws[:, ind2])
    max1 = maximum(draws[:, ind1])
    max2 = maximum(draws[:, ind2])
    max12 = max(max1, max2)
    min1 = minimum(draws[:, ind1])
    min2 = minimum(draws[:, ind2])
    min12 = min(min1, min2)
    padding = (max12 - min12) / 50.
    start = min12 - padding
    stop  = max12 + padding
    if (subspec(m) == "ss24" && occursin("full_incZLB", savename) &&
        string(DSGE.detexify(m.parameters[ind1].key)) == "psi3")
        start -= .1
    elseif (subspec(m) == "ss24" && occursin("full_incZLB", savename) &&
            string(DSGE.detexify(m.parameters[ind1].key)) == "rho")
        stop += .1
    elseif (subspec(m) == "ss24" && occursin("full_preZLB", savename) &&
            string(DSGE.detexify(m.parameters[ind1].key)) == "psi3")
        start -= .1
    elseif (subspec(m) == "ss24" && occursin("full_preZLB", savename) &&
            string(DSGE.detexify(m.parameters[ind1].key)) == "rho")
        stop += .1
    end
    @show quantile(draws[:, ind1], [0.05, 0.95])
    @show quantile(draws[:, ind2], [0.05, 0.95])

    # p = histogram(1:size(draws, 1), draws[:, ind1], color = :blue, alpha = 0.5, label = "Pre 1990",
    #               #title = title,
    #               bins = 100, normalize = :pdf)
    p = histogram(1:size(draws, 1), draws[:, ind1], color = RGB(55. / 255., 126. / 255., 184. / 255.), alpha = 0.5, label = "Pre 1990",
                  #title = title,
                  bins = 100, normalize = :pdf)
    # histogram!(p, 1:size(draws, 1), draws[:, ind2], color = :red, alpha = 0.5, label = "Post 1990", bins = 100, normalize = :pdf)
    histogram!(p, 1:size(draws, 1), draws[:, ind2], color = RGB(.8941, .1020, .1098), alpha = 0.5, label = "Post 1990", bins = 100, normalize = :pdf)
    plot!(p, range(start, length = 1000, stop = stop), pdf.(get(m.parameters[ind1].prior), range(start, length = 1000, stop = stop)), label = "Prior", color = :black, left_margin = 20px, xtickfont = font(12), ytickfont = font(12))
    # vline!(p, [mode[ind1]], color = :blue)
    # vline!(p, [mode[ind2]], color = :red)
    savefig(p, savename)
end

function add_tex(fid::IOStream, savename::String)
    @printf fid "\\includegraphics[width=.5\\textwidth]{%s} \n" savename
end

function get_param_inds(m::AbstractDSGEModel)
    tv2_inds = findall(x -> occursin("r2", string(x.key)), m.parameters)
    tv_params = map(x-> split(string(x.key), "_r2")[1], m.parameters[tv2_inds])
    inds = Vector{Vector{Int}}(undef, 0)
    for param2_ind in tv2_inds
        param = split(string(m.parameters[param2_ind].key), "_r2")[1]
        tv1_ind = findall(x -> param == string(x.key), m.parameters)
        tv2_ind = param2_ind
        push!(inds, [tv1_ind; tv2_ind])
    end
    return inds
end

fid = open("$fp/../../../docs/dsge_parameters/dsge_parameters.tex", "w")
mdds = DataFrame(model = String[], ss = String[], period = String[], split = String[], mdd = Float64[])
no_draws_yet = false
for model in ["m1002"]#, "m805", "m904", "smets_wouters_orig", "smets_wouters"]
    if model == "smets_wouters_orig"
        @printf fid "\\section{SW Orig -- Regime Switching}"
    elseif model == "smets_wouters"
        @printf fid "\\section{SW Orig -- Regime Switching}"
    else
        # @printf fid "\\section{%s -- Regime Switching}" model
    end
    for ss in ["ss22", "ss24"]#, "ss24"]#["ss21", "ss22", "ss23", "ss24", "ss25", "ss26", "ss27", "ss28", "ss29", "ss41", "ss42", "ss43", "ss44"]# ["ss45", "ss46", "ss47", "ss48"] #
        # @printf fid "\\subsection{%s}" ss
        if model == "m1002"
            m = Model1002(ss)
        elseif model == "m904"
            m = Model904(ss)
        elseif model == "m805"
            m = Model805(ss)
        elseif model == "smets_wouters_orig"
            m = SmetsWoutersOrig(ss)
        elseif model == "smets_wouters"
            m = SmetsWouters(ss)
        end
        standard_spec!(m, vint, fp; fcast_date = Date(2019, 12, 31), dsid = 10021, cdid = 1, replicate_exact = replicate_exact)
        m <= Setting(:friday, "true", true, "friday", "estimation ran on friday")
        m <= Setting(:npart, "20000", true, "npart", "number of SMC particles")
        m <= Setting(:reg, "2", true, "reg", "number of regimes")
        m <= Setting(:use_population_forecast, false)
        inds = get_param_inds(m)

        for period in ["full_incZLB", "full_preZLB"]
            if period == "full_incZLB"
                reg_splits = ["900331"]
                @printf fid "\\subsubsection{Including ZLB}"
            elseif period == "full_preZLB"
                reg_splits = ["840331", "900331"]
                @printf fid "\\subsubsection{Pre ZLB}"
            end
            m <= Setting(:period, period, true, "period", "period for this exercise (first or second)")
            for reg_split in reg_splits
                m <= Setting(:date_regime2_start_text, reg_split, true, "reg2start",
                             "The text version to be saved of when regime 2 starts")
                #try
                    mdd = marginal_data_density(m)
                    push!(mdds, [model, ss, period, reg_split, mdd])

                    draws = load_draws(m, :full)
                    mode = load_draws(m, :mode)
                    if size(draws, 2) != length(m.parameters)
                        @show size(draws, 2)
                        @error "bad $ss $period $reg_split. Delete it and re-run with correct code with correct number parameters"
                    end
                    @show "saving figures"
                    figure_save_path = "$fp/../../../save/output_data/$(spec(m))/$(subspec(m))/estimate/figures/"
                    @show figure_save_path
                    if !isdir(figure_save_path)
                        mkdir(figure_save_path)
                    end
                    @show inds
                    for ind_pair in inds
                        param_label = replace(replace(replace(DSGE.detexify(m.parameters[ind_pair[1]].tex_label),
                                              "\\" => "", ), "{" => ""), "}" => "")
                        if mod(parse(Int, split(subspec(m), "ss")[2]), 2) == 1
                            prior = "diffuse"
                        else
                            prior = "standard"
                        end
                        @show ss
                        @show "$figure_save_path$(param_label)_$(period)_$(reg_split)"
                        plot_regime_hist_comp(m, draws, mode,
                                              "$figure_save_path/$(param_label)_$(period)_$(reg_split)",
                                              ind_pair[1], ind_pair[2],
                                              "$param_label, $prior,  $period, $reg_split")
                        # add_tex(fid, "$(fp)/../../../save/output_data/$(model)/$(ss)/estimate/figures/$(param_label)_$(period)_$(reg_split).png")
                    end
                    if model == "m1002" && ss == "ss22" && reg_split == "900331"
                        mtmp = Model1002("ss10")
                        ζ_p_i_22 = m.keys[:ζ_p]
                        ζ_p_r2_i_22 = m.keys[:ζ_p_r2]
                        draws_r1 = draws[:, vcat(1:ζ_p_r2_i_22-1, ζ_p_r2_i_22+1:length(m.parameters))]
                        draws_r2 = deepcopy(draws_r1)
                        draws_r2[:, ζ_p_i_22] = draws[:, ζ_p_r2_i_22]
                        kdep, kdew = draw_κ(mtmp)

                        for i in [:κ_p, :κ_w]
                            param_label = string(DSGE.detexify(i))
                            r1draws, r2draws = if i == :κ_p
                                get_κp(mtmp, draws_r1), get_κp(mtmp, draws_r2)
                                # get_κp(m, reshape(postmode_r1, 1, DSGE.n_parameters(m))),
                                # get_κp(m, reshape(postmode_r2, 1, DSGE.n_parameters(m)))
                            else
                                get_κw(mtmp, draws_r1), get_κw(mtmp, draws_r2)
                                # get_κw(m, reshape(postmode_r1, 1, DSGE.n_parameters(m))),
                                # get_κw(m, reshape(postmode_r2, 1, DSGE.n_parameters(m)))
                            end
                            p = histogram(1:size(r1draws, 1), r1draws, color = RGB(55. / 255., 126. / 255., 184. / 255.), alpha = 0.5,
                                          label = "Pre 1990", bins = 100, normalize = :pdf)
                            histogram!(p, 1:size(r2draws, 1), r2draws, color = RGB(.8941, .1020, .1098), alpha = 0.5,
                                       label = "Post 1990", bins = 100, normalize = :pdf,
                                       left_margin = 20px, xtickfont = font(12), ytickfont = font(12),
                                       right_margin = 10px)
                            # vline!([r1postmode], label = "Mode H1", color = :blue)
                            # vline!([r2postmode], label = "Mode H2", color = :red)

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
                                stop = .1
                            end
                            if i == :κ_p
                                plot!(p, kdep.x, kdep.density, label = "Prior", color = :black,
                                      xlims = (start, stop), normalize = :pdf)
                                outer = plot(kdep.x, kdep.density, color = :black)
                            else
                                plot!(p, kdew.x, kdew.density, label = "Prior", color = :black,
                                      xlims = (start, stop), normalize = :pdf)
                            end

                            savefig(p, "$(figure_save_path)/$(param_label)_$(period)_$(reg_split).png")


                            # add_tex(fid, "$(fp)/../../../save/200129/output_data/$(model)/$(subspec(m))/estimate/figures/$(param_label)_$(ZLB).png")
                        end
                    end
                # catch err
                    # if isa(err, SystemError)
                    #     nothing
                    # else
                    #     rethrow(err)
                    # end
               # end
            end
        end
    end
end

# df_to_tex_long("$fp/../../../docs/dsge_parameters/mdds_table.tex", mdds)
close(fid)
