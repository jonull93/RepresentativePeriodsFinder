function create_plots(pf::PeriodsFinder)
    result_dir = normpath(
        joinpath(
            pf.config["base_dir"], pf.config[""]["result_dir"], "plots"
        )
    )
    !isdir(result_dir) && mkdir(result_dir)

    w = [pf.w[k] for k in sort([k for k in keys(pf.w)])]
    u = [pf.u[k] for k in sort([k for k in keys(pf.u)])]

    periods_with_weight = [i for (i,k) in enumerate(sort([k for k in keys(pf.w)])) if pf.w[k]>0]
    weights = [ww for ww in w if ww > 0]

    for ts in values(pf.time_series)

        # Original
        x = [x for x in range(1,stop=length(ts.data_norm))/length(ts.data_norm)*100.]
        y = sort(ts.data_norm, rev=true)
        plot(x, y, xlim=(0, 100), dpi=300, size = (1000, 1000/28*21), label="original")

        # Reduced
        y = ts.matrix_full_norm[periods_with_weight,:]'[:]
        x = (weights * ones(1, size(ts.matrix_full_norm)[2]))'[:] / length(ts.data_norm) * 100.
        df2 = sort(DataFrame(x=x, y=y, legend="reduced"), [:y], rev=true)
        df2[!,:x] = cumsum(df2[!,:x])
        Plots.plot!(df2.x, df2.y, label="reduced")

        xaxis!("Duration [-]", 0:10:100)
        yaxis!("Curve [-]")
        title!(ts.name)

        file_pdf = joinpath(result_dir, "$(ts.name).pdf")

        savefig(file_pdf)
    end
    #######################################################################
    # Show heatmap of sorting variable
    #######################################################################
    if any(values(pf.v) .> 0) && all(size(pf.v) .< 1000)
        v = zeros(length.([pf.periods,pf.periods])...)
        v[:, pf.rep_periods] = pf.v
        p = Plots.heatmap(
            v,
            ylabel = "Sum = Ordering of periods",
            xlabel = "Sum = Weighting of rep. period",
            title = "Diagonal = Selection of rep. period",
            aspect_ratio=:equal,
            xlim=[minimum(pf.periods),maximum(pf.periods)],
            ylim=[minimum(pf.periods),maximum(pf.periods)],
        )
        savefig(p, joinpath(result_dir, "ordering_heatmap.pdf"))
    end
end
