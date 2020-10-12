# TODO: try_get_val is actually implemnted in Base as `get` - just use that...
function try_get_val(d::AbstractDict, k, default=0.0)
    if haskey(d, k)
        return d[k]
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
