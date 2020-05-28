##################################################################################
# Author:   Hanspeter Höschle
# Date:     15/06/2017
##################################################################################
mutable struct TimeSeries
    config::Dict{Any,Any}
    time_series_type::Symbol
    weight::Float64

    name::String
    data::Array{Float64,1}                      # Real data in vector form (p*t, 1)
    data_norm::Array{Float64,1}                 # Normalised data []-1,1]
    matrix_full::Array{Float64,2}               # Real data in matrix form (p,t)
    matrix_full_norm::Array{Float64,2}          # Normalised data

    matrix_bins::Array{Float64,2}               # Relative share of days in bins (p,b); sum = 100.
    matrix_bins_cumsum_day::Array{Float64,2}    # Cumulative relative share of day in bins (p,b)
    matrix_bins_cumsum::Array{Float64,1}        # Cumulative relative share of  bins (b)

    ##################################################################################
    # Default empty constructor -> used in outer constructors in TimeSeriesFunctions.jl
    ##################################################################################
    function TimeSeries()
        return new()
    end

end
