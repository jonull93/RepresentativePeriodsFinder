############################## Copyright (C) 2019  #############################
#       The content of this file is VITO (Vlaamse Instelling voor              #
#       Technologisch Onderzoek  N.V.) proprietary.                            #
################################################################################
function get_mandatory_periods(ts::TimeSeries, dft::DaysFinderTool)

    function deletecolumn(a::Array{<:Number, 2}, del_column::Int64)
        (ncolumns, nrows) = size(a)
        remaining_rows = filter(x->!(x==del_column),collect(1:ncolumns))
        return a[remaining_rows,:]
    end
    function deletecolumn(a::Array{String, 1}, del_column::Int64)
        ncolumns = length(a)
        remaining_rows = filter(x->!(x==del_column),collect(1:ncolumns))
        return a[remaining_rows]
    end

    periods = String[]
    if isdefined(ts, :config) && "mandatory_periods" in keys(ts.config)
        for method in ts.config["mandatory_periods"]
            if method == "max"
                val, idx = findmax(ts.matrix_full)
                push!(periods, dft.periods[idx[1]])
            elseif method == "min"
                val, idx = findmin(ts.matrix_full)
                push!(periods, dft.periods[idx[1]])
            elseif (method[1] == "max") && (typeof(method[2]) == Int64)
                matrix_full = copy(ts.matrix_full)
                dftperiods = copy(dft.periods)
                periodindices = []

                for i in 1:method[2]
                    val, idx = findmax(matrix_full)
                    push!(periods, dftperiods[idx[1]])
                    push!(periodindices, idx[1])
                    dftperiods = deletecolumn(dftperiods,idx[1])
                    matrix_full = deletecolumn(matrix_full,idx[1])
                end
                println("custom max calculation executed")
            end
        end
        @debug("Mandatory periods for $(ts.name): $periods")
    end
    return unique(periods)
end
