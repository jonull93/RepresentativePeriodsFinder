"""
    create_plots(pf::PeriodsFinder, result_dir::String)

Creates plots of original and aggregated duration curves and a heatmap of the ordering variable `v` if applicable and saves them in `result_dir`.
"""
function create_plots(pf::PeriodsFinder,
        result_dir::String = get_abspath_to_result_dir(pf)
    )
    mkrootdirs(result_dir)
    rep_periods = get_set_of_representative_periods(pf)
    periods = get_set_of_periods(pf)
    weights = pf.w[rep_periods]
    ntp = get_number_of_time_steps_per_period(pf)

    for (ts_name, ts) in get_set_of_time_series(pf)

        # Original
        norm_val = get_normalised_time_series_values(pf, ts_name)
        nt = get_total_number_of_time_steps(pf)
        x = [x for x in range(1,stop=nt) / nt * 100.0]
        y = sort(norm_val[:], rev=true)
        p = Plots.plot(
            x, y, xlim=(0, 100), dpi=300, size = (800, 800/28*21), 
            label="Original"
        )

        # Reduced
        y = norm_val[rep_periods,:]'[:]
        x = (weights * ones(1, ntp))'[:] / nt * 100.0
        df = sort(DataFrame(x=x, y=y, legend="reduced"), :y, rev=true)
        df[!,:x] = cumsum(df[!,:x])

        Plots.plot!(p, df.x, df.y, label="Aggregated", title=ts_name)
        xaxis!("Duration [%]", 0:10:100)
        yaxis!("Curve [-]")

        file_pdf = joinpath(result_dir, "$(ts_name).svg")
        savefig(file_pdf)
    end

    # Heatmap of ordering variable v, if v isn't too big
    if isdefined(pf, :v) && any(values(pf.v) .> 0) && all(size(pf.v) .< 1000)
        v = zeros(length.([periods, periods])...)
        v[:, rep_periods] = pf.v
        p = Plots.heatmap(
            v,
            ylabel = "Sum = Ordering of periods",
            xlabel = "Sum = Weighting of rep. period",
            title = "Diagonal = Selection of rep. period",
            aspect_ratio=:equal,
            xlim=[minimum(periods), maximum(periods)],
            ylim=[minimum(periods), maximum(periods)],
        )
        savefig(p, joinpath(result_dir, "ordering_heatmap.svg"))
    end
end
