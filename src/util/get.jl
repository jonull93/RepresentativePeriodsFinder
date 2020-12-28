# Sets
function get_set_of_time_series_names(pf::PeriodsFinder)
    return [k for (k,v) in pf.config["time_series"] if k != "default"]::Vector{String}
end

function get_set_of_time_series(pf::PeriodsFinder)
    return [k => v for (k,v) in pf.time_series if k != "default"]::Vector{String}
end

function get_set_of_bins(pf::PeriodsFinder)
    return 1:recursive_get(pf.config, 
        "method", "optimization", "duration_curve_error", "number_bins", 40
    )::UnitRange{Int64}
end

function get_set_of_mandatory_periods(pf::PeriodsFinder)
    # TODO: fix this
    return Int64[]
end

function get_set_of_periods(pf::PeriodsFinder)
    return 1:get_number_of_periods(pf)::UnitRange{Int64}
end

function get_set_of_representative_periods(pf::PeriodsFinder)
    if isdefined(pf, :u)
        return [p for p in u if p == 1]
    else
        return nothing
    end
end

function get_set_of_ordering_errors(pf::PeriodsFinder)
    opt_general = pf.config["method"]["options"]
    if haskey(opt_general, "ordering_error")
        return [k for (k,v) in opt_general["ordering_error"]]::Vector{String}
    else
        return String[]
    end
end

# Simple parameters
function get_total_number_of_time_steps(pf::PeriodsFinder)
    opt = pf.config["method"]["options"]
    return opt["total_periods"] * opt["timesteps_per_period"]
end

function get_number_of_representative_periods(pf::PeriodsFinder)
    opt = pf.config["method"]["options"]
    return opt["representative_periods"]
end

function get_number_of_periods(pf::PeriodsFinder)
    opt = pf.config["method"]["options"]
    return opt["total_periods"]
end

function get_number_of_timesteps_per_period(pf::PeriodsFinder)
    opt = pf.config["method"]["options"]
    return opt["timesteps_per_period"]
end

function get_sampling_time(pf::PeriodsFinder)
    opt = pf.config["method"]["options"]
    return eval(Meta.parse(get(opt, "sampling_time", "Hour(1)")))
end

function get_bin_interval_values(pf::PeriodsFinder, ts_name::String)
    bins = get_set_of_bins(pf)
    x_ts = get_normalised_time_series_values(pf, ts_name)
    return bin_interval_values = range(
        minimum(x_ts), maximum(x_ts), length=length(bins)+1
    )
end

function get_error_term_weights(pf::PeriodsFinder)
    opt = pf.config["method"]["optimization"]
    opt_general = pf.config["method"]["options"]
    ord_errs = get_set_of_ordering_errors(pf)
    weights = Dict{String,Float64}()

    for err in ["time_series_error", "duration_curve_error"]
        weights[err] = recursive_get(opt, err, "weight", 0.0)
    end

    for err in ord_errs
        weights[err] = opt_general["ordering_error"][err]["weight"]
    end

    return weights
end

function get_time_series_weights(pf::PeriodsFinder)
    S = get_set_of_time_series(pf)
    return Dict(
        k => get(meta(ta), "weight", 0.0)
        for (k,ta) in S
    )
end

function get_ordering_error_function(pf::PeriodsFinder, err_name::String)
    opt_general = pf.config["method"]["options"]
    ord_err_opt = opt_general["ordering_error"]
    return errorfunc = 
        eval(Meta.parse(ord_err_opt[err_name]["function"], raise=true))
end

# Parameters

"""
    get_normalised_time_series_values(pf::PeriodsFinder)

Returns the normalised time series values from which representative periods are selected. Assigns this to `pf.x` if this has not been done already.
"""
function get_normalised_time_series_values(pf::PeriodsFinder)
    for (ts_name, ta) in get_set_of_time_series(pf)
        pf.x[ts_name] = get_normalised_time_series_values(pf, ta)
    end
    return pf.x
end

function get_normalised_time_series_values(pf::PeriodsFinder, ta::FloatTimeArray) 
    ts_name = meta(ta)["name"]
    if haskey(pf.x, ts_name) == false
        start = meta(ta)["start"]
        sampling_time = get_sampling_time(pf)
        ntt = get_total_number_of_time_steps(pf)
        timestamps = range(start, length=ntt, step=sampling_time)
        vec = values(ta[timestamps])
        
        @assert length(vec) == ntt "Less than $ntt values lie between $(timestamps[1]) and $(timestamps[end]) in $ts_name"

        vec = normalize_values(vec)

        npt = get_number_of_periods(pf)
        ntpp = get_number_of_timesteps_per_period(pf)
        return x_ts = permutedims(reshape(vec, ntpp, npt), (2,1))
    else
        return pf.x[ts_name]
    end
end

get_normalised_time_series_values(pf::PeriodsFinder, ts_name::String) = 
    get_normalised_time_series_values(pf, pf.time_series[ts_name])

"""
    get_discretised_duration_curve(pf::PeriodsFinder, ts_name::String)

Returns the duration curve of time series `ts_name` discretized into `nb` bins.
"""
function get_discretised_duration_curve(pf::PeriodsFinder, ts_name::String)
    histogram_per_period = get_histogram_per_period(pf, ts_name)
    ntt = get_total_number_of_time_steps(pf)
    L = cumsum(
        sum(histogram_per_period, dims=1) / ntt,
        dims=2
    )[:]
    return normalize_values(L)
end

function get_discretised_duration_curve(pf::PeriodsFinder)
    L = Dict(
        ts_name => get_duration_curve_parameter(pf, ts_name)
        for ts_name in get_set_of_time_series_names(pf)
    )
    return @assign L = pf.inputs
end

"""
    get_duration_curve_parameter(pf::PeriodsFinder, ts_name::String)

Returns an `np` by `nb` matrix used to approximate the aggregated duration curve of time series `ts_name` given a set of representative periods and their weights.
"""
function get_duration_curve_parameter(pf::PeriodsFinder, ts_name::String)
    histogram_per_period = get_histogram_per_period(pf, ts_name)
    np = get_number_of_periods(pf)
    A = normalize_values(
        cumsum(histogram_per_period, dims=2),
        -1/np, 1/np
    )
    return A
end

function get_duration_curve_parameter(pf::PeriodsFinder)
    A = Dict(
        ts_name => get_duration_curve_parameter(pf, ts_name)
        for ts_name in get_set_of_time_series_names(pf)
    )
    return @assign A = pf.inputs
end

"""
    get_histogram_per_period(pf::PeriodsFinder,ts_name::String)

Returns a `np`by `nb` matrix which counts the number of hours 
"""
function get_histogram_per_period(pf::PeriodsFinder,ts_name::String)
    bins = get_set_of_bins(pf)
    bin_interval_values = get_bin_interval_values(pf,ts_name)
    x_ts = get_normalised_time_series_values(pf, ts_name)
    periods = get_set_of_periods(pf)

    histogram_per_period = [
        length(
            findall(
                x -> x >= bin_interval_values[b] && x <= bin_interval_values[b+1],
                x_ts[p,:]
            )
        )
        for p in periods, b in bins
    ]

    @assert all(sum(histogram_per_period, dims=2) .== get_number_of_timesteps_per_period(pf::PeriodsFinder))

    return pf.inputs[:histogram_per_period][ts_name] = histogram_per_period
end

function get_histogram_per_period(pf::PeriodsFinder)
    histograms_per_period = Dict(
        ts_name => get_histogram_per_period(pf, ts_name)
        for ts_name in get_set_of_time_series_names(pf)
    )
    return @assign histograms_per_period = pf.inputs
end