function try_get_val(d::AbstractDict, k, default=0.0)
    if haskey(d, k)
        return d[k]
    else
        return default
    end
end
