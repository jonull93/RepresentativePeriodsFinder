function writeOutResults(dft::DaysFinderTool)
    result_dir = normpath(joinpath(dft.config["basedir"], dft.config["result_dir"]))
    if !isdir(result_dir)
        mkpath(result_dir)
    end

    df_dv = DataFrame(
                periods     = sort([k for k in dft.periods]),
                weights     = [dft.w[k] for k in sort([k for k in keys(dft.w)])],
                used_days   = [dft.u[k] for k in sort([k for k in keys(dft.u)])])
    CSV.write(joinpath(result_dir, "decision_variables.csv"), df_dv, delim=';')

    df_dv_s = deepcopy(df_dv[df_dv[!,:weights] .> 0, :])
    CSV.write(joinpath(result_dir, "decision_variables_short.csv"), df_dv_s, delim=';')

    ###########################################################################
    # Resulting time series
    ###########################################################################
    df_ts = DataFrame()
    period_idx = df_dv[!,:used_days] .== 1
    @debug("Period index of selected days: $period_idx")

    for ts in values(dft.time_series)
        df_ts[!,Symbol(ts.name)] = ts.matrix_full[period_idx,:]'[:]
    end
    CSV.write(joinpath(result_dir, "resulting_profiles.csv"), df_ts, delim=';')

    ###########################################################################
    # Save the ordering variable (if necessary)
    ###########################################################################
    order_days = try_get_val(
        dft.config, "order_days", "none"
    )
    if order_days in ["binary", "continuous"]
        IPW = [dft.v[pp,p] for p in dft.periods, pp in dft.periods]
        df = DataFrame(IPW)
        CSV.write(joinpath(result_dir, "ordering_variable.csv"), df, delim=';')
    end

    ###########################################################################
    # Copy of config-file
    ###########################################################################
    stringdata = JSON.json(dft.config)
    open(joinpath(result_dir, "config_file.json"), "w") do f
        write(f, stringdata)
    end
    return dft
end
