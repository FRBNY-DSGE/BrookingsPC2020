escape_tex(s::String) = replace(s, "_" => "\\_")
escape_tex(s::Symbol) = Symbol(escape_tex(string(s)))

function df_to_tex(filename::String, df::DataFrame)

    nrows, ncols = size(df)
    coltypes = map(colname -> eltype(df[colname]), names(df))

    open(filename, "w") do fid
        # Begin table
        @printf fid "\\begin{tabular}{"
        for i = 1:ncols
            @printf fid "c"
        end
        @printf fid "}\n"

        # Header
        for (j, colname) in enumerate(names(df))
            if j != 1
                @printf fid " & "
            end
            @printf fid "%s" escape_tex(DSGE.detexify(colname))
        end
        @printf fid " \\\\ \\hline\n"

        # Body
        for i in 1:nrows
            for j in 1:ncols
                if j != 1
                    @printf fid " & "
                end
                if coltypes[j] in [String, Symbol]
                    @printf fid "%s" "\$"
                    @printf fid "%s" escape_tex(DSGE.detexify(df[i, j]))
                    @printf fid "%s" "\$"
                else
                    @show df[i, j]
                    @printf fid "%0.3f" df[i, j]
                end
            end
            @printf fid " \\\\\n"
        end

        # End table
        @printf fid "\\end{tabular}\n"
    end

    println("Wrote $filename")
end

function df_to_tex_long(filename::String, df::DataFrame)

    nrows, ncols = size(df)
    coltypes = map(colname -> eltype(df[colname]), names(df))

    open(filename, "w") do fid
        # Begin table
        @printf fid "\\begin{longtable}{"
        for i = 1:ncols
            @printf fid "c"
        end
        @printf fid "}\n"

        # Header
        for (j, colname) in enumerate(names(df))
            if j != 1
                @printf fid " & "
            end
            @printf fid "%s" escape_tex(DSGE.detexify(colname))
        end
        @printf fid " \\\\ \\hline\n"

        # Body
        for i in 1:nrows
            for j in 1:ncols
                if j != 1
                    @printf fid " & "
                end
                if coltypes[j] in [String, Symbol]
                    @printf fid "%s" "\$"
                    @printf fid "%s" escape_tex(DSGE.detexify(df[i, j]))
                    @printf fid "%s" "\$"
                else
                    @show df[i, j]
                    @printf fid "%0.3f" df[i, j]
                end
            end
            @printf fid " \\\\\n"
        end

        # End table
        @printf fid "\\end{longtable}\n"
    end

    println("Wrote $filename")
end


function df_to_tex2(filename::String, df::DataFrame)

    nrows, ncols = size(df)
    coltypes = map(colname -> eltype(df[colname]), names(df))

    open(filename, "w") do fid
        # Begin table
        @printf fid "\\begin{tabular}{"
        for i = 1:ncols
            @printf fid "c"
        end
        @printf fid "}\n"

        # Header
        for (j, colname) in enumerate(names(df))
            if j != 1
                @printf fid " & "
            end
            @printf fid "%s" escape_tex(DSGE.detexify(colname))
        end
        @printf fid " \\\\ \\hline\n"

        # Body
        for i in 1:nrows
            for j in 1:ncols
                if j != 1
                    @printf fid " & "
                end
                if j == 1
                    @printf fid "%s" df[i,j]
                else
                    @printf fid "%0.3f" df[i, j]
                end
            end
            @printf fid " \\\\\n"
        end

        # End table
        @printf fid "\\end{tabular}\n"
    end

    println("Wrote $filename")
end
