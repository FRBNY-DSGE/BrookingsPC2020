"""
```
plot_4line(m::AbstractDSGEModel, dfs::OrderedDict{Symbol, DataFrame},
           var::Symbol, input_type::Symbol, cond_type::Symbol;
           addtl_fn_text::String = "", save_directory::String = figurespath(m, "forecast"))
```
plots the modal impulse responses corresponding to half 1 estimation (H1), half 2 estimation (H2),
and various counterfactual parametrizations.
"""
function plot_4line(m::AbstractDSGEModel, dfs::OrderedDict{Symbol, DataFrame},
                    var::Symbol, input_type::Symbol, cond_type::Symbol;
                    addtl_fn_text::String = "", save_directory::String = figurespath(m, "forecast"))
    detexify_title = typeof(Plots.backend()) == Plots.GRBackend
    if var in collect(keys(m.observables))
        the_title = DSGE.describe_series(m, var, :obs, detexify = detexify_title)
    else
        the_title = DSGE.describe_series(m, var, :pseudo, detexify = detexify_title)
    end
    p = plot(1:20, dfs[:H1][!,var], label = "Pre 1990", color = :blue, #title = the_title,
             legend = :bottomright, linewidth = 2)
    plot!(p, 1:20, dfs[:H2][!,var], label = "Post 1990", color = :red,
          linewidth = 2)
    plot!(p, 1:20, dfs[:zetap][!,var], label = "Counterfactual Slope", color = :black,
          linewidth = 2)
    plot!(p, 1:20, dfs[:MPparam][!,var], label = "Counterfactual Policy", color = :black, linestyle = :dash,
          linewidth = 2)
    if haskey(dfs, :ψπ)
        plot!(p, 1:20, dfs[:ψπ][!,var], label = "CF MP psi_pi", color = :purple)
    end
    if haskey(dfs, :Rule)
        plot!(p, 1:20, dfs[:Rule][!,var], label = "CF MP pi = pi*", color = :pink)
    end

    outfile = get_forecast_filename(save_directory, filestring_base(m),
                                    input_type, cond_type,
                                    Symbol(:line4_, DSGE.detexify(var),
                                           addtl_fn_text),
                                    forecast_string = "all",
                                    fileformat = DSGE.plot_extension())
    DSGE.save_plot(p, outfile)
end
