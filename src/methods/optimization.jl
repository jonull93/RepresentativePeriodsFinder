function make_periods_finder_model!(
        pf::PeriodsFinder, optimizer=Cbc.Optimizer
    )
    m = pf.m = JuMP.Model(optimizer)
    m.ext = Dict(
        :variables => Dict{Symbol,Any}(),
        :constraints => Dict{Symbol,Any}(),
        :expressions => Dict{Symbol,Any}(),
        :objective => nothing # Redefine later
    )
    opt = pf.config["method"]["optimization"] # optimization specifc options

    # Strict variables
    make_selection_variable!(pf, m)
    make_weighting_variable!(pf, m)
    make_ordering_variable!(pf, m)

    # Dependent error variables
    if haskey(opt, "duration_curve_error")
        make_duration_curve_error_variable!(pf, m)
        define_duration_curve_error_variable!(pf, m)
    else # else set to 0
        m.ext[:variables][:duration_curve_error] = EmptyContainer{Float64}(0.0)
    end

    if haskey(opt, "time_series_error")
        make_time_series_error_variable!(pf, m)
        define_time_series_error_variable!(pf, m)
    else # else set to 0
        m.ext[:variables][:time_series_error] = EmptyContainer{Float64}(0.0)
    end

    # Fix variables
    fix_selection_variable!(pf, m)

    # Constraints
    limit_number_of_selected_periods!(pf, m)
    restrict_non_zero_weights_to_selected_periods!(pf, m)
    enforce_yearly_equivalent_duration!(pf, m)
    restrict_non_zero_orderings_to_selected_periods!(pf, m)
    restrict_linear_combination_of_representative_periods!(pf, m)
    enforce_minimum_weight!(pf, m)

    # Objective
    formulate_objective!(pf, m)

    return m
end

function optimize_periods_finder_model!(pf::PeriodsFinder, m::JuMP.Model = pf.m)
    optimize!(pf.m)
    stat = termination_status(pf.m)

    if stat in [MOI.OPTIMAL, MOI.TIME_LIMIT] && has_values(pf.m)
        # Save selection variable as is
        pf.u = round.(var_value((pf.m.ext[:variables][:u])))

        # Save all entries of w
        rep_periods = get_set_of_representative_periods(pf)
        w = var_value(pf.m.ext[:variables][:w])
        if length(w) == length(rep_periods)
            pf.w = zeros(get_number_of_periods(pf))
            pf.w[rep_periods] = w
        elseif length(w) == get_number_of_periods(pf)
            pf.w = w
        else
            error("Length of w, $(length(w)) is unexpected.")
        end

        # Only save non-zero entries v to save space
        v = var_value(pf.m.ext[:variables][:v])
        if size(v, 2) != length(rep_periods)
            pf.v = v[:,rep_periods]
        elseif size(v, 2) == length(rep_periods)
            pf.v = v
        else
            error("Size of v, $(size(v)) is unexpected.")
        end
    else
        @show stat
    end

    # Return solution status
    return stat
end



function make_selection_variable!(pf::PeriodsFinder, m::JuMP.Model)
    rep_periods = get_set_of_representative_periods(pf)
    mandatory_periods = get_set_of_mandatory_periods(pf)
    if length(rep_periods) == length(mandatory_periods)
        return m.ext[:variables][:u] = mandatory_periods
    else
        return m.ext[:variables][:u] = @variable(m, 
            [j in rep_periods], binary=true, base_name="u"
        )
    end
end

function make_weighting_variable!(pf::PeriodsFinder, m::JuMP.Model)
    rep_periods = get_set_of_representative_periods(pf)
    opt = pf.config["method"]["optimization"]
    integer_bool = get(opt, "integral_weights", false)
    return m.ext[:variables][:w] = @variable(m, 
        [j in rep_periods], integer=integer_bool, 
        base_name="w", lower_bound=0.0
    )
end

function make_ordering_variable!(pf::PeriodsFinder, m::JuMP.Model)
    periods = get_set_of_periods(pf)
    rep_periods = get_set_of_representative_periods(pf)
    opt = pf.config["method"]["optimization"]
    bin_bool = get(opt, "binary_ordering", false)
    return m.ext[:variables][:v] = @variable(m, 
        [i in periods, j in rep_periods], 
        binary=bin_bool, base_name="v", lower_bound=0.0
    )
end

function make_duration_curve_error_variable!(pf::PeriodsFinder, m::JuMP.Model)
    S = get_set_of_time_series_names(pf)
    bins = get_set_of_bins(pf)
    return m.ext[:variables][:duration_curve_error] = @variable(m,
        [s in S, b in bins], base_name = "dc_err", lower_bound = 0.0
    )
end

function make_time_series_error_variable!(pf::PeriodsFinder, m::JuMP.Model)
    S = get_set_of_time_series_names(pf)
    rep_periods = get_set_of_representative_periods(pf)
    time_steps = get_set_of_time_steps(pf)
    return m.ext[:variables][:time_series_error] = @variable(m,
        [s in S, i in rep_periods, t in time_steps], 
        base_name = "dc_err", lower_bound = 0.0
    )
end

function fix_selection_variable!(pf::PeriodsFinder, m::JuMP.Model)
    u = m.ext[:variables][:u]
    periods = get_set_of_periods(pf)
    mandatory_periods = get_set_of_mandatory_periods(pf)
    if eltype(u) <: JuMP.VariableRef
        for p in mandatory_periods
            fix(u[p], true, force=true)
        end
    end
    return u
end

function set_selection_variable_initial_guess(pf::PeriodsFinder, m::JuMP.Model)
    # TODO: fill in
end

function limit_number_of_selected_periods!(
        pf::PeriodsFinder, m::JuMP.Model
    )
    rep_periods = get_set_of_representative_periods(pf)
    N_rep = get_number_of_representative_periods(pf)
    u = m.ext[:variables][:u]
    return m.ext[:constraints][:limit_num_periods] = 
        @constraint(m,
            sum(u[j] for j in rep_periods) <= N_rep
        )
end

function restrict_non_zero_weights_to_selected_periods!(
        pf::PeriodsFinder, m::JuMP.Model
    )
    opt = pf.config["method"]["optimization"]
    equal_weights_bool = get(opt, "equal_weights", false)
    N_total = get_number_of_periods(pf)
    N_rep = get_number_of_representative_periods(pf)
    rep_periods = get_set_of_representative_periods(pf)
    u = m.ext[:variables][:u]
    w = m.ext[:variables][:w]
    if equal_weights_bool == true
        con = [
            @constraint(m, w[j] == u[j] * N_total / N_rep)
            for j in rep_periods
        ]
    else
        con = [
            @constraint(m, w[j] <= u[j] * N_total)
            for j in rep_periods
        ]
    end
    return m.ext[:constraints][:restrict_non_zero_weights] = con
end

function enforce_yearly_equivalent_duration!(pf::PeriodsFinder, m::JuMP.Model)
    rep_periods = get_set_of_periods(pf)
    N_total = get_number_of_periods(pf)
    w = m.ext[:variables][:w]
    return m.ext[:constraints][:equivalent_yearly_duration] = @constraint(m,
        sum(w[j] for j in rep_periods) == N_total
    )
end

function restrict_non_zero_orderings_to_selected_periods!(
        pf::PeriodsFinder, m::JuMP.Model
    )
    periods = get_set_of_periods(pf)
    rep_periods = get_set_of_representative_periods(pf)
    u = m.ext[:variables][:u]
    v = m.ext[:variables][:v]
    con = [@constraint(m, v[i,j] <= u[j]) for i in periods, j in rep_periods]
    return m.ext[:constraints][:restrict_ordering] = con
end

function restrict_linear_combination_of_representative_periods!(
        pf::PeriodsFinder, m::JuMP.Model
    )
    opt = pf.config["method"]["optimization"]
    linear_comb = get(opt, "linear_combination_representative_periods", "sum_to_one")
    if linear_comb == "sum_to_one"
        periods = get_set_of_periods(pf)
        rep_periods = get_set_of_representative_periods(pf)
        v = m.ext[:variables][:v]
        con = [
            @constraint(m, sum(v[i,j] for j in rep_periods) == 1)
            for i in periods
        ]
        return m.ext[:constraints][:linear_combination_representative_periods] = con
    end
    return nothing
end

function enforce_minimum_weight!(pf::PeriodsFinder, m::JuMP.Model)
    opt = pf.config["method"]["optimization"]
    rep_periods = get_set_of_representative_periods(pf)
    min_weight = get(opt, "minimum_weight", 0.0)
    @fetch w, u = m.ext[:variables]
    return m.ext[:constraints][:minimum_weight] = [
        @constraint(m, w[j] >= u[j] * min_weight)
        for j in rep_periods
    ]
end

function define_duration_curve_error_variable!(pf::PeriodsFinder, m::JuMP.Model)
    dc_err_opt = pf.config["method"]["optimization"]["duration_curve_error"]
    type = get(dc_err_opt, "type", "squared")
    S = get_set_of_time_series_names(pf)
    bins = get_set_of_bins(pf)
    rep_periods = get_set_of_representative_periods(pf)
    A = get_duration_curve_parameter(pf)
    L = get_discretised_duration_curve(pf)
    N_total = get_number_of_periods(pf)
    @fetch duration_curve_error, w = m.ext[:variables]
    if type == "squared"
        duration_curve_error = @constraint(m,
            [s in S, b in bins],
            duration_curve_error[s][b] == 
            sum(w[j] / N_total * A[s][j,b] for j in rep_periods)
        )
    elseif type == "absolute"
        duration_curve_error_eq1 = @constraint(m, 
            [s in S, b in bins],
            duration_curve_error[s,b] >= + L[s][b] - sum(
                w[j] / N_total * A[s][j,b] for j in rep_periods
            )
        )
        duration_curve_error_eq2 = @constraint(m, 
            [s in S, b in bins],
            duration_curve_error[s,b] >= - L[s][b] + sum(
                w[j] / N_total * A[s][j,b] for j in rep_periods
            )
        )
        @assign duration_curve_error_eq1, duration_curve_error_eq2 = m.ext[:constraints]
    end
    return m
end

function define_time_series_error_variable!(pf::PeriodsFinder, m::JuMP.Model)
    return m
end

function formulate_objective!(pf::PeriodsFinder, m::JuMP.Model)
    opt = pf.config["method"]["optimization"]
    dc_err_type = recursive_get(opt, "duration_curve_error", "type", "absolute")
    ordering_errors = get_set_of_ordering_errors(pf)
    objective_term_weights = get_error_term_weights(pf)
    time_series_weights = get_time_series_weights(pf)
    periods = get_set_of_periods(pf)
    rep_periods = get_set_of_representative_periods(pf)
    bins = get_set_of_bins(pf)
    S = get_set_of_time_series_names(pf)

    ordering_error_expressions = Dict{String,Any}(
        ord_err => ordering_error_expression!(pf, m, ord_err)
        for ord_err in ordering_errors
    )

    @fetch duration_curve_error, time_series_error = m.ext[:variables]
    dc_err_func(x) = dc_err_type == "absolute" ? x^2 : x

    # @assert keys(objective_term_weights) == keys(objective_term_indices) == keys(objective_terms)

    obj = @objective(m, Min, 
        sum(
            time_series_weights[s] * (
                + objective_term_weights["duration_curve_error"] * sum(
                    dc_err_func(duration_curve_error[s,b]) for b in bins
                )
                + objective_term_weights["time_series_error"] * sum(
                    time_series_error[s,i,j] for i in periods, j in rep_periods
                )
                + sum(
                    objective_term_weights[ord_err] * sum(
                        ordering_error_expressions[ord_err][s][i,j]
                        for i in periods, j in periods
                    )
                    for ord_err in ordering_errors
                )
            )
            for s in S
        )
    )
    # for ot in keys(objective_terms), s in S
    #     for idx in Iterators.product(objective_term_indices[ot]...)
    #         try
    #             objective_term_weights[ot] * time_series_weights[s] * objective_terms[ot][s][idx...]
    #         catch e
    #             @show ot
    #             @show s
    #             @show idx
    #             @show objective_terms[ot][s][idx...]
    #             @show objective_term_weights[ot]
    #             @show time_series_weights[s]
    #             throw(e)
    #         end
    #     end
    # end

    return m.ext[:objective] = obj
end

function ordering_error_expression!(
        pf::PeriodsFinder, m::JuMP.Model, err_name::String
    )
    errorfunc = pf.inputs[:ordering_error_functions][err_name] 
    periods = get_set_of_periods(pf)
    rep_periods = get_set_of_representative_periods(pf)
    S = get_set_of_time_series_names(pf)
    x = get_normalised_time_series_values(pf)
    v = m.ext[:variables][:v]

    # @show errorfunc( [1,2], [2,3])

    ord_err_expr = Dict(
        s => [
            @expression(m, v[i,j] * errorfunc(x[s][i,:], x[s][j,:]))
            for i in periods, j in rep_periods
        ]
        for s in S
    )

    @assert all([all(length.(ord_err_expr[s]) .== 1) for s in S]) "$err_name must return scalar value."
    
    return m.ext[:expressions][Symbol(err_name)] = ord_err_expr
end

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
    dc_norm_weight = dft.N_total_periods * length(dft.time_steps) / length(dft.bins)
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
        reconstructed_time_series[c=dft.curves, i=dft.periods, t=dft.time_steps]
    )
    @constraint(m,
        [c = dft.curves, i = dft.periods, t = dft.time_steps],
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
    dc_norm_weight = dft.N_total_periods * length(dft.time_steps) / length(dft.bins)

    @objective(m, Min,
        sum(
            dft.WEIGHT_DC[c] * (
                + dc_weight * dc_norm_weight * sum(
                    (reconstructed_duration_curve[c,b] - dft.L[c,b])^2
                    for b in dft.bins
                )
                + ts_weight * sum(
                    (reconstructed_time_series[c,i,t] - dft.time_series[c].matrix_full_norm[i,t])^2
                    for i in dft.periods, t in dft.time_steps
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
        [c = dft.curves, i = P, t = dft.time_steps],
        sum(v[i,j] * dft.time_series[c].matrix_full_norm[j,t] for j = RP)
    )

    # Define error
    @variable(m, ts_error[c=dft.curves, i=P, t=dft.time_steps] >= 0)
    @constraint(m,
    	[c = dft.curves, i = P, t = dft.time_steps],
    	ts_error[c,i,t] >=
            + reconstructed_time_series[c,i,t]
    		- dft.time_series[c].matrix_full_norm[i,t]
    )
    @constraint(m,
    	[c = dft.curves, i = P,t = dft.time_steps],
    	ts_error[c,i,t] >= -(
            + reconstructed_time_series[c,i,t]
    		- dft.time_series[c].matrix_full_norm[i,t]
        )
    )

    # Define objective
    dc_weight = try_get_val(dft.config, "duration_curve_error_weight", 1.0)
    ts_weight = try_get_val(dft.config, "time_series_error_weight", 1.0)
    dc_norm_weight = dft.N_total_periods * length(dft.time_steps) / length(dft.bins)

    @objective(m, Min,
        sum(
            dft.WEIGHT_DC[c] * (
                + dc_weight * dc_norm_weight * sum(
                    duration_curve_error[c,b] for b in dft.bins
                ) + ts_weight * sum(
                    ts_error[c,i,t] for i in P, t in dft.time_steps
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
            reconstructed_time_series[c=dft.curves, i=dft.periods, t=dft.time_steps]
        )
        @constraint(m,
            [c = dft.curves, i = dft.periods, t = dft.time_steps],
            reconstructed_time_series[c,i,t]
            ==
            sum(v[i,j] * dft.time_series[c].matrix_full_norm[j,t] for j = dft.periods)
        )
    else
        reconstructed_time_series = EmptyContainer(Float64)
    end

    # Define objective
    dc_norm_weight = dft.N_total_periods * length(dft.time_steps) / length(dft.bins)

    @objective(m, Min,
        sum(
            dft.WEIGHT_DC[c] * (
                + dc_weight * dc_norm_weight * sum(
                    (reconstructed_duration_curve[c,b] - dft.L[c,b])^2
                    for b in dft.bins
                )
                + ts_weight * sum(
                    (reconstructed_time_series[c,i,t] - dft.time_series[c].matrix_full_norm[i,t])^2
                    for i in dft.periods, t in dft.time_steps
                )
            ) for c in dft.curves
        )
    )

    return m
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
                    for t in dft.time_steps))
                for p in dft.periods, pp in dft.periods
            ] for c in dft.curves
        )
    elseif type == "average"
        TSEMD = Dict(
            c => [
                errorfunc.(sum(
                    + dft.time_series[c].matrix_full_norm[p,t]
                    - dft.time_series[c].matrix_full_norm[pp,t]
                    for t in dft.time_steps
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

