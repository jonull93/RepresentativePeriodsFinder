# Sets
function get_set_of_time_series_names(pf::PeriodsFinder)
    return [k for (k,v) in pf.config["time_series"] if k != "default"]::Vector{String}
end

function get_set_of_time_series(pf::PeriodsFinder)
    S = [k => v for (k,v) in pf.time_series if k != "default"]
    return S::Vector{Pair{String,FloatTimeArray}}
end

function get_set_of_bins(pf::PeriodsFinder)
    bins = 1:recursive_get(pf.config, 
        "method", "optimization", "duration_curve_error", "number_bins", 40
    )
    return bins::UnitRange{Int64}
end

function get_set_of_mandatory_periods(pf::PeriodsFinder)
    # TODO: fix this

    mandatory_periods = [Int64(i) for i in 
        recursive_get(pf.config, "method", "options", "mandatory_periods", Int64[])
    ]
    @assert length(mandatory_periods) <= get_number_of_representative_periods(pf) 
    return mandatory_periods::Vector{Int64}
end

function get_set_of_intermediate_periods(pf::PeriodsFinder)
    intermediate_periods = [Int64(i) for i in 
        recursive_get(pf.config, "method", "clustering", 
            "intermediate_periods", Int64[]
        )
    ]
    if isempty(intermediate_periods) == false
        @assert length(intermediate_periods) < get_number_of_periods(pf)
        @assert maximum(intermediate_periods) <= get_number_of_periods(pf)
        @assert minimum(intermediate_periods) >= get_number_of_representative_periods(pf)
        @assert length(unique(intermediate_periods)) == length(intermediate_periods)
    end
    return intermediate_periods::Vector{Int64}
end

function get_set_of_periods(pf::PeriodsFinder)
    periods = 1:get_number_of_periods(pf)
    return periods::UnitRange{Int64}
end

function get_set_of_time_steps(pf::PeriodsFinder)
    periods = 1:get_number_of_time_steps_per_period(pf)
    return periods::UnitRange{Int64}
end

"""
    get_set_of_representative_periods(pf::PeriodsFinder)

If representative periods have not been defined yet, returns `get_set_of_periods`.
"""
function get_set_of_representative_periods(pf::PeriodsFinder)
    if isdefined(pf, :u) && isempty(pf.u) == false &&
            length(pf.u) == get_number_of_periods(pf)
        return [index for (index, value) in enumerate(pf.u) if value == true]
    else
        return get_set_of_periods(pf)
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

function has_ordering_error(pf::PeriodsFinder)
    return (
        isempty(get_set_of_ordering_errors(pf)) == false ||
        has_time_series_error(pf) == true
    )
end

function get_total_number_of_time_steps(pf::PeriodsFinder)
    opt = pf.config["method"]["options"]
    return opt["total_periods"] * opt["time_steps_per_period"]
end

function get_number_of_representative_periods(pf::PeriodsFinder)
    opt = pf.config["method"]["options"]
    return opt["representative_periods"]::Int64
end

function get_number_of_periods(pf::PeriodsFinder)
    opt = pf.config["method"]["options"]
    return opt["total_periods"]::Int64
end

function get_number_of_time_steps_per_period(pf::PeriodsFinder)
    opt = pf.config["method"]["options"]
    return opt["time_steps_per_period"]::Int64
end

function get_sampling_time(pf::PeriodsFinder)
    opt = pf.config["method"]["options"]
    st = get(opt, "sampling_time", Hour(1))
    typeof(st) <: AbstractString ? st = eval(Meta.parse(st)) : nothing
    return st::TimePeriod
end

function get_bin_interval_values(pf::PeriodsFinder, ts_name::String)
    bins = get_set_of_bins(pf)
    x_ts = get_normalised_time_series_values(pf, ts_name)
    return bin_interval_values = range(
        minimum(x_ts), maximum(x_ts), length=length(bins)+1
    )
end

function get_error_term_weights(pf::PeriodsFinder)
    opt_general = pf.config["method"]["options"]
    ord_errs = get_set_of_ordering_errors(pf)
    weights = Dict{String,Float64}()

    if haskey(pf.config["method"], "optimization")
        opt = pf.config["method"]["optimization"]
        for err in ["time_series_error", "duration_curve_error"]
            weights[err] = recursive_get(opt, err, "weight", 0.0)
        end
    end

    for err in ord_errs
        weights[err] = opt_general["ordering_error"][err]["weight"]
    end

    return weights
end

function has_time_series_error(pf::PeriodsFinder)
    opt = pf.config["method"]
    return haskey(opt, "optimization") && haskey(opt["optimization"], "time_series_error")
end

function has_duration_curve_error(pf::PeriodsFinder)
    opt = pf.config["method"]
    return haskey(opt, "optimization") && haskey(opt["optimization"], "duration_curve_error")
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

function get_abspath_to_result_dir(pf::PeriodsFinder)
    rel_dir = recursive_get(pf.config, "results", "result_dir", "plots")
    return abspath(joinpath(pf.config["base_dir"], rel_dir))
end

function get_selection_method(pf::PeriodsFinder)
    methods = [k for k in keys(pf.config["method"]) if k != "options"]
    @assert length(methods) == 1 "Remove either clustering or optimisation entry in .yaml config file."
    return first(methods) 
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
        ntpp = get_number_of_time_steps_per_period(pf)
        return x_ts = permutedims(reshape(vec, ntpp, npt), (2,1))
    else
        return pf.x[ts_name]
    end
end

get_normalised_time_series_values(pf::PeriodsFinder, ts_name::String) = 
    get_normalised_time_series_values(pf, pf.time_series[ts_name])

# TODO: reduce the number of lines of code
function get_time_series_values(pf::PeriodsFinder)
    vals = Dict{String,Array{Float64,2}}()
    for (ts_name, ta) in get_set_of_time_series(pf)
        start = meta(ta)["start"]
        sampling_time = get_sampling_time(pf)
        ntt = get_total_number_of_time_steps(pf)
        timestamps = range(start, length=ntt, step=sampling_time)
        vec = values(ta[timestamps])
        npt = get_number_of_periods(pf)
        ntpp = get_number_of_time_steps_per_period(pf)
        vals[ts_name] = permutedims(reshape(vec, ntpp, npt), (2,1))
    end
    return vals
end


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
        ts_name => get_discretised_duration_curve(pf, ts_name)
        for ts_name in get_set_of_time_series_names(pf)
    )
    @pack! pf.inputs = L
    return L
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
        -1, 1
    )
    return A
end

function get_duration_curve_parameter(pf::PeriodsFinder)
    A = Dict(
        ts_name => get_duration_curve_parameter(pf, ts_name)
        for ts_name in get_set_of_time_series_names(pf)
    )
    @pack! pf.inputs = A
    return A
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

    @assert all(sum(histogram_per_period, dims=2) .== get_number_of_time_steps_per_period(pf::PeriodsFinder))
    
    return recursive_set(pf.inputs, :histogram_per_period, 
        ts_name, histogram_per_period; collection_type=Dict{Union{String,Symbol},Any}
    )
end

function get_histogram_per_period(pf::PeriodsFinder)
    histograms_per_period = Dict(
        ts_name => get_histogram_per_period(pf, ts_name)
        for ts_name in get_set_of_time_series_names(pf)
    )
    @pack! pf.inputs = histograms_per_period
    return histograms_per_period
end

function get_synthetic_time_series(pf::PeriodsFinder, ts_name::String)
    return get_synthetic_time_series(pf, pf.time_series[ts_name])
end

function get_synthetic_time_series(pf::PeriodsFinder, ta::FloatTimeArray)
    start = meta(ta)["start"]
    sampling_time = get_sampling_time(pf)
    ntt = get_total_number_of_time_steps(pf)
    timestamps = range(start, length=ntt, step=sampling_time)
    vec = values(ta[timestamps])
    np = get_number_of_periods(pf)
    ntpp = get_number_of_time_steps_per_period(pf)
    array = permutedims(reshape(vec, ntpp, np), (2,1))
    rep_periods = get_set_of_representative_periods(pf)
    
    v = pf.v
    data = Float64[]
    for i in 1:np
        for t in 1:ntpp
            push!(data, v[i,:]' * array[rep_periods,t])
        end
    end
    
    return TimeArray(
        timestamp(ta[timestamps]), 
        hcat(values(ta[timestamps]), data), 
        [:original, :synthetic], meta(ta)
    )
end

# OTHER
function get_educated_guess_for_ordering_variable(pf::PeriodsFinder)
    u = pf.u
    w = pf.w
    v = fill(0.0, length(pf.u), length(pf.u))
    i = 1
    for j in 1:length(u)
        u[j] == false && continue
        n_periods = Int(round(w[j]))
        i_end = minimum([i+n_periods, length(pf.u)])
        v[i:i_end,j] .= 1
        i = i + n_periods + 1
    end
    return v
    # TODO: This "educated guess" could be improved to better deal with rounding
end