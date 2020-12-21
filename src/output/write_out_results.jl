function writeOutResults(
        dft::PeriodsFinder,
        relResultDir = dft.config["result_dir"],
    )
    result_dir = normpath(
        joinpath(dft.config["base_dir"], relResultDir)
    )
    if !isdir(result_dir)
        mkpath(result_dir)
    end

    df_dv = DataFrame(
        periods     = dft.periods,
        weights     = dft.w,
        used_days   = dft.u
    )
    CSV.write(joinpath(result_dir, "decision_variables.csv"), df_dv, delim=',')

    df_dv_s = deepcopy(sort(df_dv[df_dv[!,:used_days] .> 0, :]))
    CSV.write(joinpath(result_dir, "decision_variables_short.csv"), df_dv_s, delim=';')

    ###########################################################################
    # Resulting time series
    ###########################################################################
    df_ts = DataFrame()
    period_idx = df_dv[!,:used_days] .== 1
    # @debug("Period index of selected days: $period_idx")

    for ts in values(dft.time_series)
        df_ts[!,Symbol(ts.name)] = ts.matrix_full[period_idx,:]'[:]
    end
    CSV.write(joinpath(result_dir, "resulting_profiles.csv"), df_ts, delim=',')

    ###########################################################################
    # Save the ordering variable
    ###########################################################################
    if any(values(dft.v) .> 0)
        v = [
            dft.v[i,j] for i in dft.periods, j in 1:length(dft.rep_periods)
        ]
        permIdx = sortperm(dft.rep_periods)
        v = v[:,permIdx]
        df = DataFrame(v, Symbol.(dft.rep_periods[permIdx]))
        CSV.write(joinpath(result_dir, "ordering_variable.csv"), df, delim=',')
    end

    ###########################################################################
    # Copy of config-file
    ###########################################################################
    YAML.write_file(joinpath(result_dir, "config_file.yaml"), dft.config)
    # stringdata = JSON.json(dft.config)
    # open(joinpath(result_dir, "config_file.json"), "w") do f
    #     write(f, stringdata)
    # end

    ###########################################################################
    # Objective value and bound to yaml file
    ###########################################################################
    if isdefined(dft, :m)
        obj_val = objective_value(dft.m)
        obj_bound = objective_bound(dft.m)
        optStatus = Dict(
            "objective_value" => obj_val,
            "objective_bound" => obj_bound,
            "optimality_gap" => (obj_val - obj_bound)/obj_val * 100,
            "solution_method" => dft.config["solver"]["Method"]
        )
    else
        optStatus = Dict(
            "solution_method" => dft.config["solver"]["Method"]
        )
    end
    YAML.write_file(joinpath(result_dir, "optimisation_results.yaml"), optStatus)

    return dft
end