function time_period_clustering(dft::PeriodsFinder)
    # Notation is based on the following paper:
    # [1] S. Pineda and J. M. Morales, “Chronological time-period clustering for optimal capacity expansion planning with storage,” IEEE Trans. Power Syst., vol. 33, no. 6, pp. 7162–7170, 2018.

    type = get(dft.config["solver"], "Method", "hierarchical clustering")
    N = dft.N_total_periods # Initial number of clusters / periods
    NC = dft.N_total_periods # Current number of clusters
    NCD = dft.N_representative_periods # Desired number of cluster
    ICM = [i for i in dft.periods] # Set of periods which are the mediod of a cluster -> indices are clusters, elements are periods which are mediods of that cluster
    IP2C = [i for i in dft.periods] # A set which maps the periods to their clusters -> indices are periods, elements are clusters
    IC2P = [[i] for i in dft.periods] # A set which maps clusters to their periods -> indices are clusters, periods are elements
    x = [
        [
            Float32(dft.WEIGHT_DC[c] * dft.time_series[c].matrix_full_norm[i,t])
            for c in dft.curves, t in dft.time_steps
        ][:]
        for i in dft.periods
    ] # Cluster features vector
    xbar = copy(x) # at least initially

    # TODO: Possible performance boosts:
    # 1. Changing the order of the for loops to be j then i ->
    # be careful with this though.
    # 2. Using @inbounds -> avoids check in the array bounds when running
    # 3. Using Threads.@threads (might not speed things up though), generally bad idea
    # 5. Can also try out @simd apparently?
    # 6. Get rid of if statements in loops - just make sure i != j

    # Define dissimilarity matrix
    # Since we only care about the (unordered) pairs,
    # we only define the upper / lower triangle of the matrix
    # For hierarchical clustering, D is lower triangular and we compute 1/D
    # For CTPC, D is upper triangular with only diagonals defined
    D = fill(Inf32, NC, NC)
    # NOTE: Can also define D as UpperTriangular
    # but this leads to 0 entries which can be annoying when
    # trying to find the minimum
    Threads.@threads for j in 1:NC
        if type == "chronological time period clustering"
            istart = j
            iend = min(NC, j + 1)
        elseif type == "hierarchical clustering"
            istart = j
            iend = NC
        end
        for i in istart:iend
            if i != j
                num = Float32(2 * length(IC2P[i]) * length(IC2P[j]))
                den = Float32(length(IC2P[i]) + length(IC2P[j]))
                @inbounds D[i,j] = num / den * sum((xbar[i][:] .- xbar[j][:]).^2)
            end
        end
    end
    if type == "hierarchical clustering"
        # Do this so I can use findmax(D) and return lowest dissimilarity
        D = LowerTriangular(1 ./ D)
    end

    # Mandatory periods enforcement
    mandatoryPeriods = get(dft.config["solver"], "mandatoryPeriods", Int64[])
    replaceVal = if type == "hierarchical clustering"
        Float32(0.0)
    elseif type == "chronological time period clustering"
        Inf32
    end
    for i in mandatoryPeriods
        D[i,:] .= replaceVal
    end
    for j in mandatoryPeriods
        D[:,j] .= replaceVal
    end

    # Correct for mandatory periods
    NCD = NCD - length(mandatoryPeriods)
    @assert iszero(NCD) == false "Can't run clustering if all periods are fixed."
    # TODO: could eventually replace this with a reasonable default.
    # TODO: currently mandatory periods don't represent other periods!
    # could be annoying in the future
    ermsg = "Please allow for at least 1 period to be determined by the clustering algorithm."
    @assert NCD > 0 eval(ermsg)
    numAdj = get_number_of_clusters_of_adjacent_values(mandatoryPeriods)
    if type == "chronological time period clustering"
        ermsg =
            """
            Too many periods fixed for this method. Please remove at least
            $(NCD - numAdj + 2) periods (I think...) 
            """
        @assert NCD > numAdj + 2 eval(ermsg)
    end

    # Get intermediate periods
    intermediatePeriods = get(dft.config["solver"], "intermediatePeriods", Int64[])

    # Run clustering
    while NC - length(mandatoryPeriods) > NCD
        if mod(NC, 100) == 0
            @debug "Number of clusters: $NC"
        end

        # Find the two "closest" mediods
        # This takes up the bulk of the calculation time!
        if type == "chronological time period clustering"
            # Only need to search along a diagonal of D!
            D_vec = [
                @inbounds D[j + 1,j] for j in 1:NC - 1 # ≈ 0.1 secs
            ]
            (val, idx) = findmin(D_vec) # ≈ 10^-5 secs

            # Convert this into an index for D (instead of D_vec)
            idx = CartesianIndex(idx + 1, idx) # ≈ 0.1 secs
        else
            (val, idx) = findmax(D) # ≈ 0.5 secs
        end

        # Show if it's merging a cluster which is mandatory
        if IC2P[idx[1]] ∈ mandatoryPeriods || IC2P[idx[2]] ∈ mandatoryPeriods
            @warn "Something's probably wrong"
            @show findmax(D)
            @show NC
            @show idx[1]
            @show idx[2]
        end

        # Merge the two clusters. Do this by appending the periods of the
        # "second" cluster to that of the "first" cluster
        IC2P[idx[1]] = append!(IC2P[idx[1]], IC2P[idx[2]])

        # Now we need to define the cluster mediod
        # To do this, need the centroid of the cluster and then need to
        # calculate dissimilarity of each period w.r.t to the cluster centroid
        # The period in the cluster which minimises this dissimalarity
        # becomes the next cluster mediod

        # Recalculate the cluster centroid
        sumterm = sum(hcat(x[IC2P[idx[1]][:]]...), dims=2)[:]
        xbar[idx[1]] = sumterm / length(IC2P[idx[1]])

        # Find the merged clusters mediod as follows:
        # Calculate the dissimilarity for each element in the cluster
        DC = [sum((xbar[idx[1]] .- el).^2) for el in x[IC2P[idx[1]][:]]]

        # Find the period which minimises dissimilarity
        (val, idxClust) = findmin(DC)

        # Set this as the cluster mediod
        ICM[idx[1]] = dft.periods[IC2P[idx[1]][idxClust]]

        # Get rid of the "second" cluster in the cluster
        deleteat!(IC2P, idx[2])
        deleteat!(xbar, idx[2])
        deleteat!(ICM, idx[2])

        # Update the Periods2Cluster set
        for i = 1:length(IC2P)
            IP2C[IC2P[i][:]] .= i
        end

        # Now we update the dissimilarity matrix for the next iteration
        # Remove second cluster row and column from the dissimilarity matrix
        D_idx = [i for i in 1:NC if i != idx[2]]
        D = D[D_idx, D_idx]

        # Decrease number of clusters
        NC += -1

        # Change idx[1] in case it was the last cluster
        NC < idx[1] ? idx = CartesianIndex(NC, idx[2]) : nothing

        # Recompute the dissimilarity for the new cluster w.r.t to the other
        # Clusters
        if type == "chronological time period clustering"
            i_index = [
                i for i in max(1, idx[1] - 1):min(NC, idx[1])
            ]
            j_index = max.(1, i_index .- 1)
        elseif type == "hierarchical clustering"
            i_index = vcat(
                [idx[1] for i = 1:idx[1] - 1],
                [i for i = idx[1]:NC]
            )
            j_index = vcat(
                [j for j = 1:idx[1] - 1],
                [idx[1] for j = idx[1]:NC]
            )
        end
        @assert length(i_index) == length(j_index)
        D_idx = [
            (i_index[k], j_index[k]) for k = 1:length(i_index)
            if isempty(intersect(IC2P[i_index[k]], mandatoryPeriods))
            && isempty(intersect(IC2P[j_index[k]], mandatoryPeriods))
            && i_index[k] != j_index[k]
        ]
        for (i, j) in D_idx
            # if i <= 0 || j <= 0
            #     @show i_index
            #     @show j_index
            #     @show NC
            #     return error("Invalid indexing somehow")
            # end
            num = 2 * length(IC2P[i]) * length(IC2P[j])
            den = length(IC2P[i]) + length(IC2P[j])
            if type == "hierarchical clustering"
                @inbounds D[i,j] = Float32(
                    1 ./ (num / den * sum((xbar[i][:] .- xbar[j][:]).^2))
                )
            else
                @inbounds D[i,j] = Float32(
                    num / den * sum((xbar[i][:] .- xbar[j][:]).^2)
                )
            end
        end

        if NC in intermediatePeriods
            calculate_rep_periods_from_clusters!(dft, ICM, IC2P, IP2C)
            @warn "Currently I overwrite the saved stuff, ish sad"
            writeOutResults(
                dft, joinpath(dft.config["result_dir"], string(NC))
            )
        end
    end

    # Save the results
    calculate_rep_periods_from_clusters!(dft, ICM, IC2P, IP2C)
    return ICM, IC2P, IP2C, D
end

function calculate_rep_periods_from_clusters!(
        dft::PeriodsFinder,
        ICM,
        IC2P,
        IP2C
    )
    dft.u = [i in ICM ? 1 : 0 for i in dft.periods]
    dft.w = [i in ICM ? length(IC2P[IP2C[i]]) : 0 for i in dft.periods]
    dft.v = zeros(length(dft.periods), length(dft.periods))
    for i in dft.periods
        j = # I'm pretty sure this line is unnecessary
        dft.v[i,ICM[IP2C[i]]] = 1
    end
    dft.rep_periods = sort(ICM)
    dft.v = dft.v[:, dft.rep_periods]
    @warn "Modified dft, make sure this doesn't bugger up the clustering"
    return dft
end