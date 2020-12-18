# TODO: try_get_val is actually implemnted in Base as `get` - just use that...
function try_get_val(d::AbstractDict, k, default=0.0)
    if haskey(d, k)
        return d[k]
    else
        return default
    end
end

function config_get(d1::AbstractDict, d2::AbstractDict, k::String, default::Any)
    if haskey(d1, k)
        return d1[k]
    elseif haskey(d2, k)
        return d2[k]
    else
        return default
    end
end

function getVariableValue(x::JuMP.Containers.DenseAxisArray)
    return value.(x).data
end

function getVariableValue(x::JuMP.Containers.SparseAxisArray)
    return value.(x).data
end

function getVariableValue(x::AbstractArray)
    return x
end

function getVariableValue(x::Array{VariableRef,1})
    return value.(x)
end

function getVariableValue(x::Array{VariableRef,2})
    return value.(x)
end

function get_number_of_clusters_of_adjacent_values(x::AbstractVector)
    y = sort(x)
    numAdj = 1
    for j in 2:length(y)
        if y[j-1] + 1 == y[j]
            numAdj += 0
        else
            numAdj += 1
        end
    end
    return numAdj
end

# function get_mandatory_periods(ts::TimeSeries, dft::PeriodsFinder)
#     periods = []
#     if isdefined(ts, :config) && "mandatory_periods" in keys(ts.config)
#         for method in ts.config["mandatory_periods"]
#             if method == "max"
#                 val, idx = findmax(ts.matrix_full)
#                 push!(periods, dft.periods[idx[1]])
#             elseif method == "min"
#                 val, idx = findmin(ts.matrix_full)
#                 push!(periods, dft.periods[idx[1]])
#             end
#         end
#         @debug("Mandatory periods for $(ts.name): $periods")
#     end
#     return unique(periods)
# end
