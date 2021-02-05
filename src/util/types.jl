"""
    FloatTimeArray = TimeArray{Float64,1,DateTime,Array{Float64,1}}

The concrete TimeArray type used in PeriodsFinder.
"""
FloatTimeArray = TimeArray{Float64,1,DateTime,Array{Float64,1}}

"""
    MissingFloatTimeArray = TimeArray{T,1,DateTime,A} where T <: Union{Missing, Float64}, A <: Array{Union{Missing, Float64},1}

A concrete TimeArray type.
"""
MissingFloatTimeArray = TimeArray{T,1,DateTime,A} where T <: Union{Missing, Float64} where A <: Array{Union{Missing, Float64},1}

"""
    mutable struct EmptyContainer{T}
        x::T

Returns a specified value `x` no matter the indexing.
"""
mutable struct EmptyContainer{T}

    x::T

    EmptyContainer{T}(x) where T <: Any = new(x)

    function EmptyContainer()
        EmptyContainer(nothing)
    end

    function Base.getindex(ec::EmptyContainer, args...)
        return ec.x
    end

    function Base.print(io::IO, ec::EmptyContainer)
        print(io, "An EmptyContainer with field $(ec.x)")
    end

    function Base.show(io::IO, ec::EmptyContainer)
        print(io, ec)
    end

    function Base.length(ec::EmptyContainer)
        return length(ec.x)
    end
end
