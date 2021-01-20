"""
    PeriodsFinder(config_file::String; populate_entries::Bool=false)
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
    inputs::Dict{Any,Any} # Any saved inputs (sets, parameters, etc)

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
        println(io, "Configuration file: \n\t$(pf.config_file)")
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

    # Instantiate inputs
    ord_err_names = get_set_of_ordering_errors(pf)
    pf.inputs = Dict(
        :ordering_error_functions => Dict{String,Function}(
            ord_err => get_ordering_error_function(pf, ord_err)
            for ord_err in ord_err_names
            # Do this now to avoid "new world" issues
        )
    )

    # pf.time_series = Dict()
    # pf.curves = Array{String,1}()
    # for ts_config in pf.config["time_series"]
    #     println("-"^100)

    #     ts = TimeSeries(pf, ts_config)
    #     addTimeSeries!(pf, ts)

    #     println("$(ts.name) added")

    #     ##################################################################################
    #     # Dynamics profiles
    #     ##################################################################################
    #     if pf.config["dynamics_method"] == "all_1step"
    #         println("-"^100)

    #         ts_diff = TimeSeries(pf, ts, :all_1step)
    #         addTimeSeries!(pf, ts_diff)

    #         println("$(ts_diff.name) added")
    #     end
    # end

    # ##################################################################################
    # # Parameters
    # ##################################################################################
    # pf.WEIGHT_DC                  = Dict()
    # pf.A                          = Dict()
    # pf.L                          = Dict()

    # pf.AREA_TOTAL                 = Dict()
    # pf.AREA_TOTAL_DAY             = Dict()

    # pf.N_representative_periods   = try_get_val(pf.config, "number_days", 8)
    # pf.N_total_periods            = try_get_val(pf.config, "number_days_total", 365)

    # ##################################################################################
    # # Sets
    # ##################################################################################
    # pf.bins =     ["b" * string(b, pad=pad) for b in range(1, stop=pf.config["number_bins"])]

    # lenP = pf.N_total_periods
    # pf.periods = 1:lenP

    # # TODO: this should just be timeseries length / N_periods
    # # But this makes an assumption on the length of the timeseries, so eh...
    # lenT = try_get_val(pf.config, "time_steps_per_period", 24)
    # pf.time_steps = 1:lenT

    # ##################################################################################
    # # Correlation
    # ##################################################################################
    # if pf.config["correlation_method"] == "all"
    #     ts_basic = [ts.name for ts in values(pf.time_series) if ts.time_series_type == :basic]
    #     for combi in combinations(ts_basic, 2)
    #         println("-"^100)
    #         ts1 = pf.time_series[combi[1]]
    #         ts2 = pf.time_series[combi[2]]
    #         ts = TimeSeries(pf, ts1, ts2)
    #         addTimeSeries!(pf, ts)

    #         println("$(ts.name) added")
    #     end
    # end

    # # misc dictionary
    # pf.misc = Dict()

    return pf
end