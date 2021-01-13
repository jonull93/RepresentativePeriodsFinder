"""
    write_out_results(pf::PeriodsFinder, result_dir::String)

Writes out the representative periods selected by `pf` to `.csv` files in `result_dir`.
"""
function write_out_results(
        pf::PeriodsFinder,
        result_dir::String = get_abspath_to_result_dir(pf),
    )
    result_dir = get_abspath_to_result_dir(pf)
    mkrootdirs(result_dir)

    # Selection and weights
    df_dv = DataFrame(
        periods     = get_set_of_periods(pf),
        weights     = pf.w,
        used_days   = pf.u
    )
    CSV.write(joinpath(result_dir, "decision_variables.csv"), df_dv, delim=',')
    df_dv_s = deepcopy(sort(df_dv[df_dv[!,:used_days] .> 0, :]))
    CSV.write(joinpath(result_dir, "decision_variables_short.csv"), df_dv_s, delim=';')

    # Resulting time series
    df_ts = DataFrame()
    rep_periods = get_set_of_representative_periods(pf)
    periods = get_set_of_periods(pf)
    norm_val = get_normalised_time_series_values(pf)
    for (ts_name, ts) in get_set_of_time_series(pf)
        df_ts[!,Symbol(ts_name)] = norm_val[ts_name][rep_periods,:]'[:]
    end
    CSV.write(joinpath(result_dir, "resulting_profiles.csv"), df_ts, delim=',')

    # Save the ordering variable
    if any(values(pf.v) .> 0)
        v = [
            pf.v[i,j] for i in periods, j in 1:length(rep_periods)
        ]
        permIdx = sortperm(rep_periods)
        v = v[:,permIdx]
        df = DataFrame(v, Symbol.(rep_periods[permIdx]))
        CSV.write(joinpath(result_dir, "ordering_variable.csv"), df, delim=',')
    end

    # Copy of config-file
    YAML.write_file(joinpath(result_dir, "config_file.yaml"), pf.config)

    # Objective value and bound to yaml file
    if isdefined(pf, :m)
        obj_val = objective_value(pf.m)
        obj_bound = objective_bound(pf.m)
        optStatus = Dict(
            "objective_value" => obj_val,
            "objective_bound" => obj_bound,
            "optimality_gap" => (obj_val - obj_bound)/obj_val * 100
        )
    else
        optStatus = Dict(
            "solution_method" => pf.config["solver"]["Method"]
        )
    end
    YAML.write_file(joinpath(result_dir, "optimization_results.yaml"), optStatus)

    return pf
end

"""
    write_out_results(pf::PeriodsFinder, result_dir::String)

Writes out a synthetic time series as determined by `pf` to a `.csv` file in `result_dir`. A synthetic time series is the same length as the original time series but composed only of representative periods.
"""
function write_out_synthetic_timeseries(pf::PeriodsFinder,
        result_dir::String = get_abspath_to_result_dir(pf)
    )
    return nothing
end