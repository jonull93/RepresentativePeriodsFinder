function config_get(d1::AbstractDict, d2::AbstractDict, k::String, default::Any)
    if haskey(d1, k)
        return d1[k]
    elseif haskey(d2, k)
        return d2[k]
    else
        return default
    end
end

function mkrootdirs(dir::String)
    dirVec = split(dir, "/")
    dd = "/"
    for d in dirVec[2:end]
        dd = joinpath(dd, d)
        if isdir(dd) == false
            mkdir(dd)
        end
    end
end

function recursive_get(d, args...)
    if length(args) > 1
        if haskey(d, args[1])
            recursive_get(d[args[1]], args[2:end]...)
        else
            return args[end]
        end
    else
        return d # if length(args) > 1 then d is the value you want
    end
end

function recursive_set(d, args...)
    if length(args) > 1
        if haskey(d, args[1]) == false
            d[args[1]] = typeof(d)()
        end
        recursive_set(d[args[1]], args[2:end]...)
    end
    return last(args)
end

function normalize_values(A, lb=-1.0, ub=1.0)
    max_val = maximum(A)
    min_val = minimum(A)
    A = (A .- min_val) ./ (max_val - min_val) .* (ub - lb) .+ lb
    return replace(A, NaN => 0.0)
end

var_value(x::JuMP.Containers.DenseAxisArray) = value.(x).data
var_value(x::JuMP.Containers.SparseAxisArray) = value.(x).data
var_value(x::AbstractArray) = x
var_value(x::Array{VariableRef,1}) = value.(x)
var_value(x::Array{VariableRef,2}) = value.(x)

import Base.length
Base.length(expr::JuMP.GenericAffExpr) = 1

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

# TODO: the below can be done using the UnPack package!

"""
    @fetch x, y, ... = d

Assign mapping of `:x` and `:y` in `d` to `x` and `y` respectively.
"""
macro fetch(expr)
    (expr isa Expr && expr.head == :(=)) || error("please use @fetch with the assignment operator (=)")
    keys, dict = expr.args
    values = if keys isa Expr
        Expr(:tuple, [:($dict[$(Expr(:quote, k))]) for k in keys.args]...)
    else
        :($dict[$(Expr(:quote, keys))])
    end
    esc(Expr(:(=), keys, values))
end

"""
    @assign x, y, ... = d

Assign variables `x` and `y` to entries `:x` and `:y` in `d` respectively.
"""
macro assign(expr)
    (expr isa Expr && expr.head == :(=)) || error("please use @fetch with the assignment operator (=)")
    variables, dict = expr.args
    dictkeys = if variables isa Expr
        Expr(:tuple, [:($dict[$(Expr(:quote, k))]) for k in variables.args]...)
    else
        :($dict[$(Expr(:quote, variables))])
    end
    return esc(Expr(:(=), dictkeys, variables))
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
