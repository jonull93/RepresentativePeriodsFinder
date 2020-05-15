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
