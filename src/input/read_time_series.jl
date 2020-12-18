function add_time_series!(pf::PeriodsFinder, ts_name::String)

    ts_dict = pf.config["time_series"]

    @assert ts_name in keys(ts_dict)
    @assert ts_name != "default" "Name `default` is reserved"

    source = joinpath(
        pf.config["base_dir"],
        config_get(ts_dict[ts_name], ts_dict["default"], "source", "")
    )
    val_col = get(ts_dict[ts_name], "value_column", "")
    delim = config_get(ts_dict[ts_name], ts_dict["default"], "delim", ',')[1]
    timestamp_col = get(ts_dict[ts_name], "timestamp_column", "")
    format = config_get(ts_dict[ts_name], ts_dict["default"], "format", "")
    types = Dict(val_col => Float64, timestamp_col => DateTime)
    df = CSV.read(source, DataFrame; delim=delim, 
        types = Dict(val_col => Float64), dateformat=format
    )

    if isempty(timestamp_col)
        @debug "Assume an hourly resolution starting year 1970 for $ts_name"
        start_date = DateTime(1970, 1, 1)
        end_date = start_date + Hour(size(df,1) - 1)
        timestamps = start_date:Hour(1):end_date
        pf.time_series[ts_name] = TimeArray(timestamps, Array(df[!,val_col]))
    else
        data = eltype(df[!,val_col])[]
        timestamps = DateTime[]
        for row in eachrow(df)
            if ismissing(row[timestamp_col]) == false
                push!(data, row[val_col])
                push!(timestamps, row[timestamp_col])
            else
                continue # Replace this with an assumption on the date and time!
            end
        end
        pf.time_series[ts_name] = TimeArray(timestamps, data)
    end

    return pf.time_series[ts_name]

    # pf.time_series[ts_name] = readtimearray(source; kwargs...)

    # push!(pf.curves, ts.name)

    # weight!(pf.WEIGHT_DC, ts)
    # area_total!(pf.AREA_TOTAL, ts)
    # area_total_days!(pf.AREA_TOTAL_DAY, pf.periods, ts)

    # cum_bin_end!(pf.A, pf.periods, pf.bins, ts)
    # cum_bin_total!(pf.L, self.bins, ts)
end