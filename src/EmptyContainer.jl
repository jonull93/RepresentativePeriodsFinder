"""
        mutable struct EmptyContainer

Returns zero(T) or nothing if T <: Nothing.
"""
mutable struct EmptyContainer

    returnType::Type # Determines what kind of zero value you return

    """
        EmptyContainer(T::Type)

    Create an EmptyContainer.
    """
    function EmptyContainer(T::Type)
        ec = new()
        if (T <: Nothing) == false
            try
                zero(T)
            catch
                error("Cannot return zero of type $T")
            end
        end
        ec.returnType = T
        return ec
    end

    function EmptyContainer()
        EmptyContainer(Nothing)
    end

    function Base.getindex(ec::EmptyContainer, args...)
        if ec.returnType <: Nothing
            return nothing
        else
            return zero(ec.returnType)
        end
    end

    function Base.print(io::IO, ec::EmptyContainer)
        print(io, "An EmptyContainer")
    end

    function Base.show(io::IO, ec::EmptyContainer)
        print(io, ec)
    end

    function Base.length(ec::EmptyContainer)
        return 0
    end
end
