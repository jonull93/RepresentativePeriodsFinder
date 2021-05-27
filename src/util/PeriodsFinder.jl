"""
    PeriodsFinder(config_file::String; populate_entries::Bool=false)

Type used for finding representative periods in a time series (see `find_representative_periods`). If `populate_entries=true` then time series are loaded into `pf`.
"""
mutable struct PeriodsFinder
    ###########################################################################
    # Configuration
    ###########################################################################
    config_file::AbstractString
    config::Dict{Any,Any}

    ###########################################################################
    # Inputs
    ###########################################################################
    time_series::Dict{String,FloatTimeArray}
    x::Dict{String,Array{Float64,2}} # Normalised time series in matrix format
    inputs::Dict{Symbol,Any} # Any saved inputs (sets, parameters, etc)
    # TODO: Change this to Dict{Union{Symbol,String},Any} - much nicer.

    ###########################################################################
    # optimization model
    ###########################################################################
    m::JuMP.Model

    ###########################################################################
    # Results
    ###########################################################################
    u::Array{Bool,1}    # Selection variable
    w::Array{Float64,1} # Weight variable (can also be integer)
    v::Array{T,2} where T <: Union{Bool,Float64} # Ordering variable (not always defined, can be float)
    # TODO: redefine v so it's Union{Bool,Float64} to potentially save memort

    function PeriodsFinder()
        pf = new()
        pf.config = Dict{Any,Any}()
        pf.time_series = Dict{String,TimeArray}()
        return pf
    end

    function PeriodsFinder(config_file::String;     
            populate_entries::Bool=false
        )
        pf = new()
        pf.config_file = config_file
        @assert isfile(config_file)
        pf.config = YAML.load(open(config_file))
        !(haskey(pf.config, "base_dir")) && (pf.config["base_dir"] = dirname(config_file))
        pf.time_series = Dict{String,TimeArray}()

        if populate_entries == true
            pf = populate_entries!(pf)
        end

        return pf
    end

    # Pretty prints of long term planning model
    function Base.print(io::IO, pf::PeriodsFinder)
        println(io, "Representative Periods Finder")
        isdefined(pf, :config_file) && println(io, "Configuration file: \n\t$(pf.config_file)")
    end

    function Base.show(io::IO, pf::PeriodsFinder)
        print(io, pf)
    end
end

function populate_entries!(pf::PeriodsFinder)
    S = get_set_of_time_series_names(pf)
    
    for ts_name in S
        @info "Adding $ts_name..."
        ta = read_time_series(pf, ts_name)
        if haskey(meta(ta), "interpolation_type")
            interpolate_missing_values!(pf, ta)
        end
        if haskey(meta(ta), "resample") && meta(ta)["resample"] == true
            resample!(pf, ta)
        end
        pf.time_series[ts_name] = TimeArray(
            DateTime.(timestamp(ta)), Float64.(values(ta)), 
            colnames(ta), meta(ta)
        )
        # TODO: 
        # 1) Create correlation time series (and weights)
        # 2) Create ramping time series (and weights)
    end
       
    pf.x = Dict{String,Array{Float64,2}}()

    # Instantiate inputs, load error functions
    ord_err_names = get_set_of_ordering_errors(pf)
    pf.inputs = Dict(
        :ordering_error_functions => Dict{String,Function}(
            ord_err => get_ordering_error_function(pf, ord_err)
            for ord_err in ord_err_names
            # Do this now to avoid "new world" issues
        )
    )

    return pf
end