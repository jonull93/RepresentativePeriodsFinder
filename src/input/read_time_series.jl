"""
    read_time_series(pf::PeriodsFinder, ts_name::String)

Reads a time series based on data in `pf.config[ts_name]`.
"""
function read_time_series(pf::PeriodsFinder, ts_name::String)

    ts_dict = pf.config["time_series"]

    @assert ts_name in keys(ts_dict)
    @assert ts_name != "default" "Name `default` is reserved"

    # Read in .csv file
    source = joinpath(
        pf.config["base_dir"],
        config_get(ts_dict[ts_name], ts_dict["default"], "source", "")
    )
    val_col = get(ts_dict[ts_name], "value_column", "")
    timestamp_col = get(ts_dict[ts_name], "timestamp_column", "")
    csv_options = Dict()
    for name in ("default", ts_name)
        if haskey(ts_dict[name], "csv_options") && typeof(ts_dict[name]["csv_options"]) <: Dict
            merge!(csv_options, ts_dict[name]["csv_options"])
        end
    end
    csv_options = Dict{Symbol,Any}(
        Symbol(k) => get(csv_options, k, "") 
        for (k,v) in csv_options if isempty(v) == false
    )
    csv_options[:types] = Dict(val_col => Float64, timestamp_col => DateTime)
    csv_options = namedtuple(collect(csv_options))
    df = CSV.read(source, DataFrame; csv_options...)

    # Checks
    for col in filter(!isempty, [val_col, timestamp_col])
        standard_err_string = """

        Options passed to CSV.read: $csv_options

        Head of data frame: $(first(df, 5))

        Tail of data frame: $(last(df, 5))
        """

        @assert col in names(df) """ 
        "$col" not a column name in $source. 
        $standard_err_string
        """

        col_type = eltype(df[!,col])
        req_type = Union{Missing,csv_options[:types][col]}
        @assert col_type <: req_type """
        "$col" is of type $(col_type) and not $(req_type). 
        $standard_err_string
        """
    end

    # Data for TimeArray type
    if isempty(timestamp_col)
        @debug "Assume an hourly resolution starting year 1970 for $ts_name"
        start_date = DateTime(1970, 1, 1)
        end_date = start_date + Hour(size(df, 1) - 1)
        timestamps = start_date:Hour(1):end_date
        data = Array(df[!,val_col])
    else
        data = eltype(df[!,val_col])[]
        timestamps = DateTime[]
        for row in eachrow(df)
            if ismissing(row[timestamp_col]) == false
                push!(data, row[val_col])
                push!(timestamps, row[timestamp_col])
            else
                continue # Replace this with an assumption on the date and time?
            end
        end
    end

    metaDict = Dict(
        k => v for (k, v) in merge(ts_dict[ts_name], ts_dict["default"])
    )
    metaDict["name"] = ts_name
    if haskey(metaDict, "start")
        start = metaDict["start"]
        if typeof(start) <: Int
            start = timestamps[start]
        else
            start = DateTime(start, DateFormat(csv_options[:dateformat]))
        end
    else
        start = timestamps[1]
    end
    metaDict["start"] = start

    return ta = TimeArray(timestamps, data, [:value], metaDict)

    # pf.time_series[ts_name] = readtimearray(source; kwargs...)

    # push!(pf.curves, ts.name)

    # weight!(pf.WEIGHT_DC, ts)
    # area_total!(pf.AREA_TOTAL, ts)
    # area_total_days!(pf.AREA_TOTAL_DAY, pf.periods, ts)

    # cum_bin_end!(pf.A, pf.periods, pf.bins, ts)
    # cum_bin_total!(pf.L, self.bins, ts)
end