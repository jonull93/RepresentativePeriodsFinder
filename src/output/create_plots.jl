"""
    create_plots(pf::PeriodsFinder, result_dir::String)

Creates plots of original and aggregated duration curves and a heatmap of the ordering variable `v` if applicable and saves them in `result_dir`.
"""
function create_plots(pf::PeriodsFinder,
        result_dir::String = get_abspath_to_result_dir(pf)
    )
    mkrootdirs(result_dir)
    rep_periods = get_set_of_representative_periods(pf)
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
        x = transpose(weights * ones(1, ntp))[:] / nt * 100.0
        df = sort(DataFrame(x=x, y=y, legend="reduced"), :y, rev=true)
        df[!,:x] = cumsum(df[!,:x])

        Plots.plot!(p, df.x, df.y, label="Aggregated", title=ts_name)
        xaxis!("Duration [%]", 0:10:100)
        yaxis!("Curve [-]")

        file_svg = joinpath(result_dir, "$(ts_name)_duration_curve.svg")
        savefig(file_svg)
    end

    # Heatmap of ordering variable v, if v isn't too big
    if isdefined(pf, :v) && any(values(pf.v) .> 0) && all(size(pf.v) .< 1000)
        p = create_ordering_heatmap(pf)
        savefig(p, joinpath(result_dir, "ordering_heatmap.svg"))
    end
    return nothing
end

"""
    create_ordering_heatmap(pf)
    
Returns a heatmap showing the mapping of representative periods to non-representative periods. 
"""
function create_ordering_heatmap(pf::PeriodsFinder)
    rep_periods = get_set_of_representative_periods(pf)
    periods = get_set_of_periods(pf)
    v = zeros(length.([periods, periods])...)
    v[:, rep_periods] = pf.v
    return Plots.heatmap(
        v,
        ylabel = "Sum = Ordering of periods",
        xlabel = "Sum = Weighting of rep. period",
        title = "Diagonal = Selection of rep. period",
        aspect_ratio=:equal,
        xlim=[minimum(periods), maximum(periods)],
        ylim=[minimum(periods), maximum(periods)],
    )
end

"""
    create_synthetic_time_series_plots(pf::PeriodsFinder,
        result_dir::String = get_abspath_to_result_dir(pf);
        timestamps
    )

Creates plots of original and synthetic time series for the ranges in `timestamp` and saves them in `result_dir`.

# Example
```julia
timestamps = Dict("Load" => DateTime(1970,1,1):Hour(1):DateTime(1970,1,2))
# result_dir not specified, since it is specified in pf.config 
create_synthetic_time_series_plots(pf; timestamps=timestamps)
```
"""
function create_synthetic_time_series_plots(pf::PeriodsFinder,
        result_dir::String = get_abspath_to_result_dir(pf);
        timestamps = Dict{String,Vector{DateTime}}()
    )
    mkrootdirs(result_dir)
    for (ts_name, ts) in get_set_of_time_series(pf)
        p = create_synthetic_time_series_plot(pf, ts_name; 
            timestamps=timestamps
        )
        savefig(p, joinpath(result_dir, "$(ts_name)_synthetic_time_series.svg"))
    end
    return nothing
end

"""
    create_synthetic_time_series_plot(pf::PeriodsFinder, ts_name; timestamps)
    
# Example
```julia
timestamps = Dict("Load" => DateTime(1970,1,1):Hour(1):DateTime(1970,1,2))
p = create_synthetic_time_series_plot(pf, "Load"; timestamps=timestamps)
display(p)
```
"""
function create_synthetic_time_series_plot(pf::PeriodsFinder,
        ts_name::String;
        timestamps = Dict{String,Vector{DateTime}}()
    )
    ts = pf.time_series[ts_name]
    tstamp = haskey(timestamps, ts_name) ? timestamps[ts_name] : timestamp(ts)
    return Plots.plot(
        get_synthetic_time_series(pf, ts)[tstamp],
        title=ts_name,
    )
end