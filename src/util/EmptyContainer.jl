"""
    mutable struct EmptyContainer{T}
        x::T

Returns a specified value `x` no matter the indexing.
"""
mutable struct EmptyContainer{T}

    x::T # Determines what kind of zero value you return

    function EmptyContainer()
        EmptyContainer(nothing)
    end

    function Base.getindex(ec::EmptyContainer, args...)
        if ec.returnType <: Nothing
            return nothing
        else
            return zero(ec.returnType)
        end
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
