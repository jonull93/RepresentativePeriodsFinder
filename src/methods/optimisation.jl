function makeDCErrorOnlyPeriodsFinderModel(
    dft::PeriodsFinder, optimizer_factory
    )
    # Model
    m = Model(optimizer_factory)
    dft.m = m

    # Variables
    @variable(m, u[j in dft.periods], Bin)
    if try_get_val(dft.config, "integral_weights", false)
        @variable(m, w[j in dft.periods], Int)
    else
        @variable(m, w[j in dft.periods] >= 0)
    end
    v = spzeros(length.([dft.periods,dft.periods])...)
    dft.misc[:u] = u
    dft.misc[:v] = v
    dft.misc[:w] = w

    # Defining error as absolute value
    @variable(m, duration_curve_error[c in dft.curves, b in dft.bins] >= 0)
    @constraint(m, duration_curve_error_eq1[c in dft.curves, b in dft.bins],
        duration_curve_error[c,b] >= + dft.L[c,b] - sum(w[j] / dft.N_total_periods * dft.A[c,j,b] for j in dft.periods)
    )
    @constraint(m, duration_curve_error_eq2[c in dft.curves, b in dft.bins],
        duration_curve_error[c,b] >= - dft.L[c,b] + sum(w[j] / dft.N_total_periods * dft.A[c,j,b] for j in dft.periods)
    )

    # User defined number of representative periods
    @constraint(m, number_periods_eq,
        sum(u[j] for j in dft.periods) == dft.N_representative_periods
    )

    # Restrict non-zero weights to selected periods
    equal_weights = try_get_val(dft.config, "equal_weights", false)
    if equal_weights == true
        @constraint(m, single_weight_eq[j in dft.periods],
            w[j] == u[j] * dft.N_total_periods / dft.N_representative_periods
        )
    else
        @constraint(m, single_weight_eq[j in dft.periods],
            w[j] <= u[j] * dft.N_total_periods
        )
    end

    # Guarantee equivalent yearly duration
    @constraint(m, total_weight_eq,
        sum(w[j] for j in dft.periods) == dft.N_total_periods
    )

    # Minimum weight
    @debug("minimum weight set to 0.1")
    @constraint(m, minimum_weight[p in dft.periods],
        w[p] >= u[p] * 0.1
    )

    # Objective
    dc_weight = try_get_val(dft.config, "duration_curve_error_weight", 1.0)
    dc_norm_weight = dft.N_total_periods * length(dft.timesteps) / length(dft.bins)
    obj = @expression(m,
        sum(
            dft.WEIGHT_DC[c] * (
                + dc_weight * dc_norm_weight * sum(
                    duration_curve_error[c,b] for b in dft.bins
                )
            ) for c in dft.curves
        )
    )
    @objective(m, Min, obj)
    return m
end

function makeSquaredErrorPeriodsFinderModel(
    dft::PeriodsFinder, optimizer_factory
    )
    # Model
    m = Model(optimizer_factory)
    dft.m = m

    # Variables
    @variable(m, u[j in dft.periods], Bin)
    if try_get_val(dft.config, "integral_weights", false)
        @variable(m, w[j in dft.periods], Int)
    else
        @variable(m, w[j in dft.periods] >= 0)
    end
    @variable(m, v[j in dft.periods, i in dft.periods] >= 0)
    dft.misc[:u] = u
    dft.misc[:v] = v
    dft.misc[:w] = w

    # Define duration curve
    @variable(m,
        reconstructed_duration_curve[c in dft.curves, b in dft.bins] >= 0
    )
    @constraint(m, duration_curve_con[c in dft.curves, b in dft.bins],
        reconstructed_duration_curve[c,b]
        == sum(w[j] / dft.N_total_periods * dft.A[c,j,b] for j in dft.periods)
    )

    # Define the time series
    @variable(m,
        reconstructed_time_series[c=dft.curves, i=dft.periods, t=dft.timesteps]
    )
    @constraint(m,
        [c = dft.curves, i = dft.periods, t = dft.timesteps],
        reconstructed_time_series[c,i,t]
        ==
        sum(v[i,j] * dft.time_series[c].matrix_full_norm[j,t] for j = dft.periods)
    )

    # Enforce v < u
    for i in dft.periods, j in dft.periods
        @constraint(m,
            v[i,j] <= u[j]
        )
    end

    # User defined number of representative periods
    @constraint(m, number_periods_eq,
        sum(u[j] for j in dft.periods) == dft.N_representative_periods
    )

    # Restrict non-zero weights to selected periods
    equal_weights = try_get_val(dft.config, "equal_weights", false)
    if equal_weights == true
        @constraint(m, single_weight_eq[j in dft.periods],
            w[j] == u[j] * dft.N_total_periods / dft.N_representative_periods
        )
    else
        @constraint(m, single_weight_eq[j in dft.periods],
            w[j] <= u[j] * dft.N_total_periods
        )
    end

    # Guarantee equivalent yearly duration
    @constraint(m, total_weight_eq,
        sum(w[j] for j in dft.periods) == dft.N_total_periods
    )

    # Minimum weight
    @debug("minimum weight set to 0.1")
    @constraint(m, minimum_weight[p in dft.periods],
        w[p] >= u[p] * 0.1
    )

    # Objective
    # Define objective
    dc_weight = try_get_val(dft.config, "duration_curve_error_weight", 1.0)
    ts_weight = try_get_val(dft.config, "time_series_error_weight", 1.0)
    dc_norm_weight = dft.N_total_periods * length(dft.timesteps) / length(dft.bins)

    @objective(m, Min,
        sum(
            dft.WEIGHT_DC[c] * (
                + dc_weight * dc_norm_weight * sum(
                    (reconstructed_duration_curve[c,b] - dft.L[c,b])^2
                    for b in dft.bins
                )
                + ts_weight * sum(
                    (reconstructed_time_series[c,i,t] - dft.time_series[c].matrix_full_norm[i,t])^2
                    for i in dft.periods, t in dft.timesteps
                )
            ) for c in dft.curves
        )
    )

    return m
end

function makePeriodsFinderModel(dft::PeriodsFinder, optimizer_factory)
    # Model
    m = Model(optimizer_factory)
    dft.m = m

    # Create variables
    v = Array{VariableRef,2}(
        undef, length(dft.periods), length(dft.periods)
    )
    u = Array{VariableRef,1}(undef, length(dft.periods))
    w = Array{VariableRef,1}(undef, length(dft.periods))

    @debug "Making ordering variable..."
    for i in dft.periods, j in dft.periods
        v[i,j] = @variable(m, binary = true)
        if mod(i - 1, 1000) == 0 && j == 1
            @debug "i = $i"
        end
    end
    for i in dft.periods
        u[i] = @variable(m, binary = true)
    end
    if try_get_val(dft.config, "integral_weights", false)
        for i in dft.periods
            w[i] = @variable(m, lower_bound = 0)
        end
    else
        for i in dft.periods
            w[i] = @variable(m, integer = true, lower_bound = 0)
        end
    end

    # Save these to PeriodsFinder
    dft.misc[:u] = u
    dft.misc[:v] = v
    dft.misc[:w] = w

    # # Hot start v
    # v_start_rand = getRandomHotStartForOrderingVariable(dft)
    # v_start = try_get_val(
    #     dft.config, "hot_start_values", v_start_rand
    # )
    # for i in dft.periods, j in dft.periods
    #     set_start_value(v[i,j], v_start[i,j])
    # end

    # If v[i,j] = 1, then day i is represented by day j
    # The sum over j (columns) is the weighting of a particular day
    # The sum over i (rows) must equal 1 - a day is represented by only one other day

    # Each day i has to be represented by another day j
    @debug "Enforcing representation of each period by a representative period..."
    for i in dft.periods
        @constraint(m,
            sum(v[i,j] for j in dft.periods) == 1
        )
    end

    # If day j is selected as a representative day then u[j] is 1
    @debug "Enforcing selection of representative period..."
    for i in dft.periods, j in dft.periods
        @constraint(m,
            v[i,j] <= u[j]
        )
    end

    linear_comb = try_get_val(
        dft.config["solver"], "LinearCombination", "sum_to_one"
    )
    if linear_comb == "sum_to_one"
        # TODO: this is actually exactly the same as the first constraint...
        # Can probably drop it!
        @constraint(m, [i = dft.periods], sum(v[i,j] for j = dft.periods) == 1)
    elseif linear_comb == "relaxed"
        nothing
    end

    # Choose only N_representative periods
    @debug "Constraint number of representative periods..."
    @constraint(m,
        sum(u[j] for j in dft.periods) == dft.N_representative_periods
    )

    # Define weightings - either in relation to the chronology variable v or seperately
    @debug "Define weightings..."
    disj_weight = try_get_val(dft.config["solver"], "DecoupledWeights", false)
    if disj_weight == true
        for j in dft.periods
            @constraint(m,
                w[j] <= u[j] * dft.N_total_periods
            )
        end
    else
        for j in dft.periods
            @constraint(m,
                w[j] == sum(v[i,j] for i in dft.periods)
            )
        end
    end

    # Equal weightings
    equal_weights = try_get_val(dft.config, "equal_weights", false)
    if equal_weights == true
        for j in dft.periods
            @constraint(m,
                w[j] == u[j] * dft.N_total_periods / dft.N_representative_periods
            )
        end
    end

    # Sum of weightings of representative periods must equal N_total_periods
    @constraint(m,
        sum(w[j] for j in dft.periods) == dft.N_total_periods
    )

    # Other constraints
    # TODO Impose particular day as a solution

    # Define objective
    @debug "Getting time series error matrix..."
    TSEMD = getTimeSeriesErrorMatrix(dft)
    @debug "Getting duration curve error matrix..."
    DCEMD = getDurationCurveErrorMatrix(dft)
    dc_weight = try_get_val(dft.config, "duration_curve_error_weight", 1.0)
    ts_weight = try_get_val(dft.config, "time_series_error_weight", 1.0)

    @debug "Making objective..."
    obj = @expression(m,
        sum(
            dft.WEIGHT_DC[c] * (
                + dc_weight * sum(
                    v[i,j] * DCEMD[c][i,j] for i in dft.periods, j in dft.periods
                )
                + ts_weight * sum(
                    v[i,j] * TSEMD[c][i,j] for i in dft.periods, j in dft.periods
                )
            ) for c in dft.curves
        )
    )
    @objective(m, Min, obj)

    # Constrain lower bound of objective
    # TODO: probably not relevant
    obj_lower_bound = try_get_val(
        dft.config, "objective_lower_bound", 0.0
    )
    @constraint(m, obj >= obj_lower_bound)

    return m
end

function makeReOrderingPeriodsFinder(
    dft::PeriodsFinder,
    optimizer_factory
    )
    println("-"^80)
    println("Re ordering days")
    println("-"^80)

    # Make model
    m = Model(optimizer_factory)
    dft.m = m

    # Get sets
    if isdefined(dft, :rep_periods)
        RP = dft.rep_periods # rep. periods j
    else
        error("Representative periods not defined, cannot reorder")
    end
    P = dft.periods # all periods i

    # Make ordering variable and weighting expression
    if try_get_val(dft.config, "integral_ordering", false)
        @variable(m, v[i=P,j=RP], Bin)
    else
        @variable(m, v[i=P,j=RP] >= 0)
    end

    if try_get_val(dft.config, "integral_weights", false)
        @variable(m, w[j in dft.periods], Int)
    else
        @variable(m, w[j in dft.periods] >= 0)
    end
    @constraint(m, [j in RP],
        w[j] == sum(v[i,j] for i in dft.periods)
    )
    dft.misc[:v] = v
    dft.misc[:w] = w
    dft.u = [i in RP ? 1 : 0 for i in dft.periods]
    dft.misc[:u] = dft.u # Indicates that selection has been done

    # Enforce same weightings if applicable
    enforce_same_weightings = try_get_val(
        dft.config, "enforce_same_weightings", false
    )
    if enforce_same_weightings == true
        @constraint(m, [j = RP],
            w[j] == dft.w[j]
        )
    end

    # Enforce equal weightings
    equal_weights = try_get_val(dft.config, "equal_weights", false)
    if equal_weights == true
        @constraint(m, single_weight_eq[j in RP],
            w[j] == u[j] * dft.N_total_periods / dft.N_representative_periods
        )
    end

    # Enforce that representative period is selected for it's own period
    @constraint(m, [j = RP], v[j,j] == 1)

    # Ensure that every day is assigned a linear combination of rep periods
    linear_comb = try_get_val(
        dft.config["solver"], "LinearCombination", "sum_to_one"
    )
    if linear_comb == "sum_to_one"
        @constraint(m, [i = P], sum(v[i,j] for j = RP) == 1)
    elseif linear_comb == "relaxed"
        nothing
    end

    # Ensure same yearly duration
    @constraint(m, sum(v[i,j] for i = P, j = RP) == dft.N_total_periods)

    # Define DC error
    @variable(m, duration_curve_error[c in dft.curves, b in dft.bins] >= 0)
    @constraint(m, error_eq1[c in dft.curves, b in dft.bins],
        duration_curve_error[c,b] >= + dft.L[c,b] - sum(w[j] / dft.N_total_periods * dft.A[c,j,b] for j in RP)
    )
    @constraint(m, error_eq2[c in dft.curves, b in dft.bins],
        duration_curve_error[c,b] >= - dft.L[c,b] + sum(w[j] / dft.N_total_periods * dft.A[c,j,b] for j in RP)
    )

    # Make expression for the reconstructed time series
    reconstructed_time_series = @expression(m,
        [c = dft.curves, i = P, t = dft.timesteps],
        sum(v[i,j] * dft.time_series[c].matrix_full_norm[j,t] for j = RP)
    )

    # Define error
    @variable(m, ts_error[c=dft.curves, i=P, t=dft.timesteps] >= 0)
    @constraint(m,
    	[c = dft.curves, i = P, t = dft.timesteps],
    	ts_error[c,i,t] >=
            + reconstructed_time_series[c,i,t]
    		- dft.time_series[c].matrix_full_norm[i,t]
    )
    @constraint(m,
    	[c = dft.curves, i = P,t = dft.timesteps],
    	ts_error[c,i,t] >= -(
            + reconstructed_time_series[c,i,t]
    		- dft.time_series[c].matrix_full_norm[i,t]
        )
    )

    # Define objective
    dc_weight = try_get_val(dft.config, "duration_curve_error_weight", 1.0)
    ts_weight = try_get_val(dft.config, "time_series_error_weight", 1.0)
    dc_norm_weight = dft.N_total_periods * length(dft.timesteps) / length(dft.bins)

    @objective(m, Min,
        sum(
            dft.WEIGHT_DC[c] * (
                + dc_weight * dc_norm_weight * sum(
                    duration_curve_error[c,b] for b in dft.bins
                ) + ts_weight * sum(
                    ts_error[c,i,t] for i in P, t in dft.timesteps
                )
            ) for c in dft.curves
        )
    )

    return m
end

function makeSquaredErrorReOrderingPeriodsFinder(
    dft::PeriodsFinder,
    optimizer_factory
    )
    println("-"^80)
    println("Re ordering days")
    println("-"^80)

    # Make model
    m = Model(optimizer_factory)
    dft.m = m

    # Get sets
    if isdefined(dft, :rep_periods)
        RP = dft.rep_periods # rep. periods j
    else
        error("Representative periods not defined, cannot reorder")
    end
    P = dft.periods # all periods i

    # Make ordering and weighting variable
    if try_get_val(dft.config, "integral_ordering", false)
        @variable(m, v[i=P,j=RP], Bin)
    else
        @variable(m, v[i=P,j=RP] >= 0)
    end

    if try_get_val(dft.config, "integral_weights", false)
        @variable(m, w[j in dft.periods], Int)
    else
        @variable(m, w[j in dft.periods] >= 0)
    end
    @constraint(m, [j in RP],
        w[j] == sum(v[i,j] for i in dft.periods)
    )
    dft.misc[:v] = v
    dft.misc[:w] = w
    dft.u = [i in RP ? 1 : 0 for i in dft.periods]
    dft.misc[:u] = dft.u # Indicates that selection has been done

    # Enforce same weightings if applicable
    enforce_same_weightings = try_get_val(
        dft.config, "enforce_same_weightings", false
    )
    if enforce_same_weightings == true
        @constraint(m, [j = RP],
            w[j] == dft.w[j]
        )
    end

    # Enforce equal weightings
    equal_weights = try_get_val(dft.config, "equal_weights", false)
    if equal_weights == true
        @constraint(m, single_weight_eq[j in RP],
            w[j] == u[j] * dft.N_total_periods / dft.N_representative_periods
        )
    end

    # Enforce that representative period is selected for it's own period
    @constraint(m, [j = RP], v[j,j] == 1)

    # Ensure that every day is assigned a linear combination of rep periods
    linear_comb = try_get_val(
        dft.config["solver"], "LinearCombination", "sum_to_one"
    )
    if linear_comb == "sum_to_one"
        @constraint(m, [i = P], sum(v[i,j] for j = RP) == 1)
    elseif linear_comb == "relaxed"
        nothing
    end

    # Ensure same yearly duration
    @constraint(m, sum(v[i,j] for i = P, j = RP) == dft.N_total_periods)

    # Get DC and TS weights
    dc_weight = try_get_val(dft.config, "duration_curve_error_weight", 1.0)
    ts_weight = try_get_val(dft.config, "time_series_error_weight", 1.0)

    # Define DC curve
    if dc_weight > 0
        @variable(m,
            reconstructed_duration_curve[c in dft.curves, b in dft.bins] >= 0
        )
        @constraint(m, duration_curve_con[c in dft.curves, b in dft.bins],
            reconstructed_duration_curve[c,b]
            == sum(w[j] / dft.N_total_periods * dft.A[c,j,b] for j in dft.periods)
        )
    else
        reconstructed_duration_curve = EmptyContainer(Float64)
    end

    # Define the time series
    if ts_weight > 0
        @variable(m,
            reconstructed_time_series[c=dft.curves, i=dft.periods, t=dft.timesteps]
        )
        @constraint(m,
            [c = dft.curves, i = dft.periods, t = dft.timesteps],
            reconstructed_time_series[c,i,t]
            ==
            sum(v[i,j] * dft.time_series[c].matrix_full_norm[j,t] for j = dft.periods)
        )
    else
        reconstructed_time_series = EmptyContainer(Float64)
    end

    # Define objective
    dc_norm_weight = dft.N_total_periods * length(dft.timesteps) / length(dft.bins)

    @objective(m, Min,
        sum(
            dft.WEIGHT_DC[c] * (
                + dc_weight * dc_norm_weight * sum(
                    (reconstructed_duration_curve[c,b] - dft.L[c,b])^2
                    for b in dft.bins
                )
                + ts_weight * sum(
                    (reconstructed_time_series[c,i,t] - dft.time_series[c].matrix_full_norm[i,t])^2
                    for i in dft.periods, t in dft.timesteps
                )
            ) for c in dft.curves
        )
    )

    return m
end

function solveGrowAlgorithmForPeriodSelection(
    dft::PeriodsFinder,
    optimizer_factory
    )
    RP = Array{Int64,1}()
    N_representative_periods = dft.N_representative_periods
    # dft.config["solver"]["Method"] = "reorder"
    for j = 1:N_representative_periods
        periodsIterator = [i for i in dft.periods if i ∉ RP]
        objVals = Array{Float64,1}(undef, length(periodsIterator))
        for k in 1:length(periodsIterator)
            dft.rep_periods = sort(vcat(RP, periodsIterator[k]))
            dft.N_representative_periods = j
            makeReOrderingPeriodsFinder(dft, optimizer_factory)
            stat = optimizePeriodsFinder(dft)
            objVals[k] = objective_value(dft.m)
        end
        sortIdx = sortperm(objVals)
        RP = sort(vcat(RP, periodsIterator[sortIdx[1]]))
    end
    return RP
end

function optimizePeriodsFinder(dft::PeriodsFinder)
    optimize!(dft.m)
    stat = termination_status(dft.m)

    # Write results to dft if optimal
    if stat in [MOI.OPTIMAL, MOI.TIME_LIMIT] && has_values(dft.m)
        # Save selection variable as is
        dft.u = round.(getVariableValue(dft.misc[:u]))

        # Save representative periods and their mapping
        dft.rep_periods = [
            idx for (idx, val) in enumerate(dft.u) if val > 0
        ]

        # NOTE: length(dft.rep_periods) != dft.N_representative_periods necessarily! Could be that optimiser chose less days than it was allowed to.

        # Transform w variable to be correct length if needs be
        w = getVariableValue(dft.misc[:w])
        if length(w) == length(dft.rep_periods)
            dft.w = zeros(dft.N_total_periods)
            dft.w[dft.rep_periods] = w
        else
            dft.w = w
        end

        # Only save the non-zero columns of v
        v = getVariableValue(dft.misc[:v])
        if size(v, 2) != length(dft.rep_periods)
            dft.v = v[:,dft.rep_periods]
        else
            dft.v = v
        end
    else
        @show stat
    end

    # Return solution status
    return stat
end

function getTimeSeriesErrorMatrix(dft::PeriodsFinder)
    type = try_get_val(
        dft.config, "time_series_error_matrix_type", "average"
    )
    squareBool = try_get_val(
        dft.config, "use_squared_error_matrix", false
    )
    errorfunc = squareBool ? x -> x^2 : x -> abs(x)
    if type == "absolute"
        TSEMD = Dict(
            c => [
                sum(errorfunc.(+ dft.time_series[c].matrix_full_norm[p,t]
                    - dft.time_series[c].matrix_full_norm[pp,t]
                    for t in dft.timesteps))
                for p in dft.periods, pp in dft.periods
            ] for c in dft.curves
        )
    elseif type == "average"
        TSEMD = Dict(
            c => [
                errorfunc.(sum(
                    + dft.time_series[c].matrix_full_norm[p,t]
                    - dft.time_series[c].matrix_full_norm[pp,t]
                    for t in dft.timesteps
                ))
                for p in dft.periods, pp in dft.periods
            ] for c in dft.curves
        )
    end
    return TSEMD
end

function getDurationCurveErrorMatrix(dft::PeriodsFinder)
    # TODO: This is a symmetric matrix, exploit that.
    squareBool = try_get_val(
        dft.config, "use_squared_error_matrix", false
    )
    errorfunc = squareBool ? x -> x^2 : x -> abs(x)
    DCEMD = Dict(
        c => [
            sum(
                errorfunc.(.+ sort(dft.time_series[c].matrix_full_norm[p,:])
                    .- sort(dft.time_series[c].matrix_full_norm[pp,:]))
            )
            for p in dft.periods, pp in dft.periods
        ] for c in dft.curves
    )

    # Periodic duration curve minimising the square error
    # DCEMD = Dict(
    #     c => [
    #         sum(
    #             abs.(
    #                 .+ sort(dft.time_series[c].matrix_full_norm[p,:])
    #                 .- sort(dft.time_series[c].matrix_full_norm[pp,:])
    #             ).^2
    #         )
    #         for p in dft.periods, pp in dft.periods
    #     ] for c in dft.curves
    # )

    # The below gives the error in the range - is a stupid idea though.
    # DCEMD = Dict(
    #     c => [
    #             abs.(
    #                 (
    #                     + maximum(dft.time_series[c].matrix_full_norm[p,:])
    #                     - minimum(dft.time_series[c].matrix_full_norm[p,:])
    #                 )
    #                 -
    #                 (
    #                     + maximum(dft.time_series[c].matrix_full_norm[pp,:])
    #                     - minimum(dft.time_series[c].matrix_full_norm[pp,:])
    #                 )
    #             )
    #         for p in dft.periods, pp in dft.periods
    #     ] for c in dft.curves
    # )
end

function getRandomHotStartForOrderingVariable(dft::PeriodsFinder)
    v_start_rand = spzeros(length.([dft.periods,dft.periods])...)
    idx = rand(dft.periods, dft.N_representative_periods - 1)
    v_start_rand[idx,idx] .= 1
    return v_start_rand
end

function fix_periods!(dft::PeriodsFinder)
    u_fix = try_get_val(dft.config, "fixed_periods", [])
    u = dft.misc[:u]
    length(u_fix) > 0 ? @info("Fixing periods $u_fix") : nothing
    for i in 1:length(u_fix)
        ui = u[u_fix[i]]
        if typeof(ui) <: JuMP.VariableRef
            fix(ui, 1, force=true)
        end
    end
end