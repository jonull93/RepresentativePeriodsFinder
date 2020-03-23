############################## Copyright (C) 2019  #############################
#       The content of this file is VITO (Vlaamse Instelling voor              #
#       Technologisch Onderzoek  N.V.) proprietary.                            #
################################################################################
function create_plots(dft::DaysFinderTool)
    result_dir = normpath(joinpath(dft.config["basedir"], dft.config["result_dir"], "plots"))
    if !isdir(result_dir)
        mkdir(result_dir)
    end

    w = [dft.w[k] for k in sort([k for k in keys(dft.w)])]
    u = [dft.u[k] for k in sort([k for k in keys(dft.u)])]

    periods_with_weight = [i for (i,k) in enumerate(sort([k for k in keys(dft.w)])) if dft.w[k]>0]
    weights = [ww for ww in w if ww > 0]

    for ts in values(dft.time_series)

        # Original
        x = [x for x in range(1,stop=length(ts.data))/length(ts.data)*100.]
        y = sort(ts.data, rev=true)
        plot(x, y, xlim=(0, 100), dpi=300, size = (1000, 1000/28*21), label="original")

        # Reduced
        y = ts.matrix_full[periods_with_weight,:]'[:]
        x = (weights * ones(1, size(ts.matrix_full)[2]))'[:] / length(ts.data) * 100.
        df2 = sort(DataFrame(x=x, y=y, legend="reduced"), [:y], rev=true)
        df2[!,:x] = cumsum(df2[!,:x])
        plot!(df2.x, df2.y, label="reduced")

        xaxis!("Duration [-]", 0:10:100)
        yaxis!("Curve [-]")
        title!(ts.name)

        file_pdf = joinpath(result_dir, "$(ts.name).pdf")

        savefig(file_pdf)
    end
    #######################################################################
    # Show heatmap of sorting variable
    #######################################################################
    v = Array{Float64,2}([dft.v[pp,p] for p in dft.periods, pp in dft.periods])
    p = heatmap(
        v,
        xlabel = "Sum = Ordering of periods",
        ylabel = "Sum = Weighting of rep. period",
        title = "Diagonal = Selection of rep. period",
        aspect_ratio=:equal,
        xlim=[minimum(dft.periods),maximum(dft.periods)],
        ylim=[minimum(dft.periods),maximum(dft.periods)],
    )
    savefig(p, joinpath(result_dir, "ordering_heatmap.pdf"))
end
