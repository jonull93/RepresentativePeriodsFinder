# function add_all_time_series!(pf)
#     for ts_name in diff(config[""], 

function add_time_series!(pf::PeriodsFinder, ts_name::String)

    ts_dict = pf.config["time_series"]

    @assert ts_name in keys(ts_dict)
    @assert ts_name != "default" "Name `default` is reserved"

    source = joinpath(
        pf.config["base_dir"],
        config_get(ts_dict[ts_name], ts_dict["default"], "source", "")
    )
    column = get(ts_dict[ts_name], "value_column", "")
    delim = config_get(ts_dict[ts_name], ts_dict["default"], "delim", ',')[1]

    df = CSV.read(source, DataFrame; delim=delim)

    timestamp = get(ts_dict[ts_name], "timestamp_column", "")

    if isempty(timestamp)
        # Assume an hourly resolution starting year 1970 for one year
        start_date = DateTime(1970, 1, 1)
        end_date = start_date + Hour(size(df,1) - 1)
        timestamps = start_date:Hour(1):end_date
        pf.time_series[ts_name] = TimeArray(timestamps, df[!,column])
    else
        format = config_get(ts_dict[ts_name], ts_dict["default"], "format", "")
        timestamps = DateTime.(df[!,timestamp], format)
        pf.time_series[ts_name] = TimeArray(timestamps, df[!,column])
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