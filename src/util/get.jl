# Sets
function get_set_of_time_series_names(pf::PeriodsFinder)
    return [k for (k,v) in pf.config["time_series"] if k != "default"]
end

function get_set_of_time_series(pf::PeriodsFinder)
    return [k => v for (k,v) in pf.time_series if k != "default"]
end

function get_set_of_bins(pf::PeriodsFinder)
    return recursive_get(pf.config, 
        "method", "optimisation", "duration_curve_error", "number_bins", 40
    )
end

# Simple parameters
function get_total_number_of_time_steps(pf::PeriodsFinder)
    return pf.config["total_periods"] * pf.config["timesteps_per_period"]
end

function get_sampling_time(pf::PeriodsFinder)
    return eval(Meta.parse(get(pf.config, "sampling_time", "Hour(1)")))
end

# Parameters
function get_normalised_time_series_values(pf::PeriodsFinder)
    if isdefined(pf, :x) == false || isempty(pf.x) == true
        for (ts_name, ta) in get_set_of_time_series(pf)
            ts_name = meta(ta)["name"]
            start = meta(ta)["start"]

            sampling_time = get_sampling_time(pf)
            ntt = get_total_number_of_time_steps(pf)
            timestamps = range(start, length=ntt, step=sampling_time)
            vec = values(ta[timestamps])
            @assert length(vec) == ntt "Less than $ntt values lie between $(timestamps[1]) and $(timestamps[end]) in $ts_name"

            vec = normalize_values(vec)

            npt = pf.config["total_periods"]
            ntpp = pf.config["timesteps_per_period"]
            pf.x[ts_name] = permutedims(reshape(vec, ntpp, npt), (2,1))
        end
    end
    return pf.x
end

function get_duration_curve_parameter(pf::PeriodsFinder)

end