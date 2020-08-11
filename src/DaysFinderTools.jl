###############################################################################
# Author:   Hanspeter Höschle
# Date:     15/06/2017
##################################################################################
mutable struct DaysFinderTool
    ##################################################################################
    # Configuration
    ##################################################################################
    config_file::AbstractString
    config::Dict{Any,Any}

    ##################################################################################
    # Sets
    ##################################################################################
    bins::Array                   # set of bins
    curves::Array                 # set of duration curves
    periods::Array{Int64,1}            # set of potential representative periods
    rep_periods::Array{Int64,1}   # set of chosen representative periods - list of integers
    timesteps::UnitRange # Set of timesteps per potential period

    ##################################################################################
    # Parameters
    ##################################################################################
    N_representative_periods::Int           # Number of representative periods to select
    N_total_periods::Int                    # Total number of repetitions required to scale up the duration of a single representative period to one year

    WEIGHT_DC::Dict                         # Relative importance of approximating duration curve c
    A::Dict            # Share of the time of day d during which the lowest value of the range corresponding to bin b of duration curve c is exceeded
    L::Dict                   # Share of the time during which the values of a time series with corresponding duration curve c exceed the lowest value of the range corresponding to bin b
    AREA_TOTAL::Dict                        # Total area under the DC of time series c
    AREA_TOTAL_DAY::Dict      # Total area under time series with DC c per day
    AREA_ERROR_TOLERANCE::Dict              # tolerance parameter for area constraint

    ###########################################################################
    # Model
    ###########################################################################
    m::JuMP.Model

    ##################################################################################
    # Time Series
    ##################################################################################
    time_series::Dict{String,TimeSeries}

    ##################################################################################
    # Results
    ##################################################################################
    u::Array{Int64,1}
    w::Array{Float64,1}
    v::Array{Float64,2}

    ###########################################################################
    # Misc dictionary
    ###########################################################################
    misc::Dict

    function DaysFinderTool(config_file::String; populate_entries::Bool = false)
        self = new()
        self.config_file = config_file
        self.config = YAML.load(open(config_file))
        self.config["basedir"] = dirname(config_file)

        if populate_entries == true
            self = populateDaysFinderTool!(self)
        end

        return self
    end

    # Pretty prints of long term planning model
    function Base.print(io::IO, dft::DaysFinderTool)
        println(io, "Representative Days Finder")
        println(io, "Configuration file: \n\t$(dft.config_file)")
    end

    function Base.show(io::IO, dft::DaysFinderTool)
        print(io, dft)
    end
end

function populateDaysFinderTool!(self::DaysFinderTool)
    pad = 3
    ##################################################################################
    # Parameters
    ##################################################################################
    self.WEIGHT_DC                  = Dict()
    self.A                          = Dict()
    self.L                          = Dict()

    self.AREA_TOTAL                 = Dict()
    self.AREA_TOTAL_DAY             = Dict()

    self.N_representative_periods   = try_get_val(self.config, "number_days", 8)
    self.N_total_periods            = try_get_val(self.config, "number_days_total", 365)

    ##################################################################################
    # Sets
    ##################################################################################
    self.bins =     ["b"*string(b, pad=pad) for b in range(1,stop=self.config["number_bins"])]

    lenP = self.N_total_periods
    self.periods = 1:lenP

    # TODO: this should just be timeseries length / N_periods
    # But this makes an assumption on the length of the timeseries, so eh...
    lenT = try_get_val(self.config, "timesteps_per_period", 24)
    self.timesteps = 1:lenT

    ##################################################################################
    # Initialize the TimeSeries
    ##################################################################################
    self.time_series = Dict()
    self.curves = Array{String,1}()
    for ts_config in self.config["time_series"]
        println("-"^100)

        ts = TimeSeries(self, ts_config)
        addTimeSeries!(self, ts)

        println("$(ts.name) added")

        ##################################################################################
        # Dynamics profiles
        ##################################################################################
        if self.config["dynamics_method"] == "all_1step"
            println("-"^100)

            ts_diff = TimeSeries(self, ts, :all_1step)
            addTimeSeries!(self, ts_diff)

            println("$(ts_diff.name) added")
        end
    end

    ##################################################################################
    # Correlation
    ##################################################################################
    if self.config["correlation_method"] == "all"
        ts_basic = [ts.name for ts in values(self.time_series) if ts.time_series_type ==:basic]
        for combi in combinations(ts_basic,2)
            println("-"^100)
            ts1 = self.time_series[combi[1]]
            ts2 = self.time_series[combi[2]]
            ts = TimeSeries(self, ts1, ts2)
            addTimeSeries!(self, ts)

            println("$(ts.name) added")
        end
    end

    # misc dictionary
    self.misc = Dict()

    return self
end

function addTimeSeries!(self::DaysFinderTool, ts::TimeSeries)
    self.time_series[ts.name] = ts
    push!(self.curves, ts.name)

    weight!(self.WEIGHT_DC, ts)
    area_total!(self.AREA_TOTAL, ts)
    area_total_days!(self.AREA_TOTAL_DAY, self.periods, ts)

    cum_bin_end!(self.A, self.periods, self.bins, ts)
    cum_bin_total!(self.L, self.bins, ts)
end

function makeDCErrorOnlyDaysFinderToolModel(
    dft::DaysFinderTool, optimizer_factory
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
        duration_curve_error[c,b] >= + dft.L[c,b] - sum( w[j] / dft.N_total_periods * dft.A[c,j,b] for j in dft.periods)
    )
    @constraint(m, duration_curve_error_eq2[c in dft.curves, b in dft.bins],
        duration_curve_error[c,b] >= - dft.L[c,b] + sum( w[j] / dft.N_total_periods * dft.A[c,j,b] for j in dft.periods)
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
    dc_norm_weight = dft.N_total_periods*length(dft.timesteps)/length(dft.bins)
    obj = @expression(m,
        sum(
            dft.WEIGHT_DC[c] *(
                + dc_weight*dc_norm_weight*sum(
                    duration_curve_error[c,b] for b in dft.bins
                )
            ) for c in dft.curves
        )
    )
    @objective(m, Min, obj)
    return m
end

function makeSquaredErrorDaysFinderToolModel(
    dft::DaysFinderTool, optimizer_factory
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
        [c=dft.curves, i=dft.periods, t=dft.timesteps],
        reconstructed_time_series[c,i,t]
        ==
        sum(v[i,j]*dft.time_series[c].matrix_full_norm[j,t] for j=dft.periods)
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
    dc_norm_weight = dft.N_total_periods*length(dft.timesteps)/length(dft.bins)

    @objective(m, Min,
        sum(
            dft.WEIGHT_DC[c] *(
                + dc_weight*dc_norm_weight*sum(
                    (reconstructed_duration_curve[c,b] - dft.L[c,b])^2
                    for b in dft.bins
                )
                + ts_weight*sum(
                    (reconstructed_time_series[c,i,t] - dft.time_series[c].matrix_full_norm[i,t])^2
                    for i in dft.periods, t in dft.timesteps
                )
            ) for c in dft.curves
        )
    )

    return m
end

function makeDaysFinderToolModel(dft::DaysFinderTool, optimizer_factory)
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
        v[i,j] = @variable(m, binary=true)
        if mod(i-1, 1000) == 0 && j == 1
            @debug "i = $i"
        end
    end
    for i in dft.periods
        u[i] = @variable(m, binary=true)
    end
    if try_get_val(dft.config, "integral_weights", false)
        for i in dft.periods
            w[i] = @variable(m, lower_bound=0)
        end
    else
        for i in dft.periods
            w[i] = @variable(m, integer=true, lower_bound=0)
        end
    end

    # Save these to DaysFinderTool
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
        @constraint(m, [i=dft.periods], sum(v[i,j] for j=dft.periods) == 1)
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
            dft.WEIGHT_DC[c] *(
                + dc_weight*sum(
                    v[i,j]*DCEMD[c][i,j] for i in dft.periods, j in dft.periods
                )
                + ts_weight*sum(
                    v[i,j]*TSEMD[c][i,j] for i in dft.periods, j in dft.periods
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

function makeReOrderingDaysFinderTool(
    dft::DaysFinderTool,
    optimizer_factory
    )
    println("-" ^ 80)
    println("Re ordering days")
    println("-" ^ 80)

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
        @constraint(m, [j=RP],
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
    @constraint(m, [j=RP], v[j,j] == 1)

    # Ensure that every day is assigned a linear combination of rep periods
    linear_comb = try_get_val(
        dft.config["solver"], "LinearCombination", "sum_to_one"
    )
    if linear_comb == "sum_to_one"
        @constraint(m, [i=P], sum(v[i,j] for j=RP) == 1)
    elseif linear_comb == "relaxed"
        nothing
    end

    # Ensure same yearly duration
    @constraint(m, sum(v[i,j] for i=P, j=RP) == dft.N_total_periods)

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
        [c=dft.curves, i=P, t=dft.timesteps],
        sum(v[i,j]*dft.time_series[c].matrix_full_norm[j,t] for j=RP)
    )

    # Define error
    @variable(m, ts_error[c=dft.curves, i=P, t=dft.timesteps] >= 0)
    @constraint(m,
    	[c=dft.curves, i=P, t=dft.timesteps],
    	ts_error[c,i,t] >=
            + reconstructed_time_series[c,i,t]
    		- dft.time_series[c].matrix_full_norm[i,t]
    )
    @constraint(m,
    	[c=dft.curves, i=P,t=dft.timesteps],
    	ts_error[c,i,t] >= -(
            + reconstructed_time_series[c,i,t]
    		- dft.time_series[c].matrix_full_norm[i,t]
        )
    )

    # Define objective
    dc_weight = try_get_val(dft.config, "duration_curve_error_weight", 1.0)
    ts_weight = try_get_val(dft.config, "time_series_error_weight", 1.0)
    dc_norm_weight = dft.N_total_periods*length(dft.timesteps)/length(dft.bins)

    @objective(m, Min,
        sum(
            dft.WEIGHT_DC[c] *(
                + dc_weight*dc_norm_weight*sum(
                    duration_curve_error[c,b] for b in dft.bins
                ) + ts_weight*sum(
                    ts_error[c,i,t] for i in P, t in dft.timesteps
                )
            ) for c in dft.curves
        )
    )

    return m
end

function makeSquaredErrorReOrderingDaysFinderTool(
    dft::DaysFinderTool,
    optimizer_factory
    )
    println("-" ^ 80)
    println("Re ordering days")
    println("-" ^ 80)

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
        @constraint(m, [j=RP],
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
    @constraint(m, [j=RP], v[j,j] == 1)

    # Ensure that every day is assigned a linear combination of rep periods
    linear_comb = try_get_val(
        dft.config["solver"], "LinearCombination", "sum_to_one"
    )
    if linear_comb == "sum_to_one"
        @constraint(m, [i=P], sum(v[i,j] for j=RP) == 1)
    elseif linear_comb == "relaxed"
        nothing
    end

    # Ensure same yearly duration
    @constraint(m, sum(v[i,j] for i=P, j=RP) == dft.N_total_periods)

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
            [c=dft.curves, i=dft.periods, t=dft.timesteps],
            reconstructed_time_series[c,i,t]
            ==
            sum(v[i,j]*dft.time_series[c].matrix_full_norm[j,t] for j=dft.periods)
        )
    else
        reconstructed_time_series = EmptyContainer(Float64)
    end

    # Define objective
    dc_norm_weight = dft.N_total_periods*length(dft.timesteps)/length(dft.bins)

    @objective(m, Min,
        sum(
            dft.WEIGHT_DC[c] *(
                + dc_weight*dc_norm_weight*sum(
                    (reconstructed_duration_curve[c,b] - dft.L[c,b])^2
                    for b in dft.bins
                )
                + ts_weight*sum(
                    (reconstructed_time_series[c,i,t] - dft.time_series[c].matrix_full_norm[i,t])^2
                    for i in dft.periods, t in dft.timesteps
                )
            ) for c in dft.curves
        )
    )

    return m
end

function solveGrowAlgorithmForPeriodSelection(
    dft::DaysFinderTool,
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
            makeReOrderingDaysFinderTool(dft, optimizer_factory)
            stat = optimizeDaysFinderTool(dft)
            objVals[k] = objective_value(dft.m)
        end
        sortIdx = sortperm(objVals)
        RP = sort(vcat(RP, periodsIterator[sortIdx[1]]))
    end
    return RP
end

function timePeriodClustering(dft::DaysFinderTool)
    # Notation is based on the following paper:
    # [1] S. Pineda and J. M. Morales, “Chronological time-period clustering for optimal capacity expansion planning with storage,” IEEE Trans. Power Syst., vol. 33, no. 6, pp. 7162–7170, 2018.

    # DEbugging stuff, ignore or delete
    # config_file = "/home/u0128861/Desktop/RepDays4Storage/rep-days-4-storage/input_data/Elia_2017.yaml"
    # dft = DaysFinderTool(config_file)
    # populateDaysFinderTool!(dft)

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
            for c in dft.curves, t in dft.timesteps
        ][:]
        for i in dft.periods
    ] # Cluster features vector
    xbar = copy(x) # at least initially

    # Define dissimilarity matrix
    # Since we only care about the (unordered) pairs, we only define the upper triangle of the matrixs
    D = fill(Inf, NC, NC)
    # NOTE: Can also define D as UpperTriangular
    # but this leads to 0 entries which is annoying when
    # trying to find the minimum
    Threads.@threads for j in 1:NC
        if type == "chronological time period clustering"
            istart = j
            iend = min(NC, j+1)
        elseif type == "hierarchical clustering"
            istart = j
            iend = NC
        end
        for i in istart:iend
            if i != j
                num = Float32(2*length(IC2P[i])*length(IC2P[j]))
                den = Float32(length(IC2P[i]) + length(IC2P[j]))
                @inbounds D[i,j] = num/den * sum((xbar[i][:] .- xbar[j][:]) .^ 2)
            end
        end
    end
    if type == "hierarchical clustering"
        # Do this so I can use findmax(D) and return lowest dissimilarity
        D = LowerTriangular(1 ./ D)
    end

    # TODO: Other possible performance boosts:
    # 1. Changing the order of the for loops to be j then i ->
    # Would need to change to LowerTriangular
    # 2. Using @inbounds -> avoids check in the array bounds when running
    # 3. Using Threads.@threads (might not speed things up though)
    # 4. Using Float32 instead of Float64
    # 5. Can also try out @simd apparently?

    intermediatePeriods = get(dft.config["solver"], "intermediatePeriods", [])

    while NC > NCD
        if mod(NC, 100) == 0
            @debug "Number of clusters: $NC"
        end

        # Find the two "closest" mediods
        # This takes up the bulk of the calculation time!
        if type == "chronological time period clustering"
            # TODO: Check this! Probably wrong
            # Only need to search along a diagonal of D!
            D_vec = [
                @inbounds D[j+1,j] for j in 1:NC-1
            ]
            (val, idx) = findmin(D_vec)

            # Convert this into an index for D (instead of D_vec)
            idx = CartesianIndex(idx + 1, idx)
        else
            (val, idx) = findmax(D)
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
        xbar[idx[1]] = sumterm/length(IC2P[idx[1]])

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
        NC < idx[1] ? idx = CartesianIndex(NC,idx[2]) : nothing

        # Recompute the dissimilarity for the new cluster w.r.t to the other
        # Clusters
        if type == "chronological time period clustering"
            i_index = [i for i in max(1, idx[1] - 1):min(NC,idx[1])]
            j_index = max.(1, i_index .- 1)
            # Double check that
        elseif type == "hierarchical clustering"
            i_index = vcat(
                [idx[1] for i = 1:idx[1]-1],
                [i for i = idx[1]:NC]
            )
            j_index = vcat(
                [j for j = 1:idx[1]-1],
                [idx[1] for j = idx[1]:NC]
            )
        end
        @assert length(i_index) == length(j_index)
        D_idx = [(i_index[k],j_index[k]) for k = 1:length(i_index)]
        for (i,j) in D_idx
            if i != j
                if i <= 0 || j <= 0
                    @show i_index
                    @show j_index
                    @show NC
                    return error("Invalid indexing somehow")
                end
                num = 2*length(IC2P[i])*length(IC2P[j])
                den = length(IC2P[i]) + length(IC2P[j])
                if type == "hierarchical clustering"
                    @inbounds D[i,j] = Float32(
                        1 ./ (num/den * sum((xbar[i][:] .- xbar[j][:]) .^ 2))
                    )
                else
                    @inbounds D[i,j] = Float32(
                        num/den * sum((xbar[i][:] .- xbar[j][:]) .^ 2)
                    )
                end
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
    return ICM, IC2P, IP2C
end

function calculate_rep_periods_from_clusters!(
        dft::DaysFinderTool,
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

function optimizeDaysFinderTool(dft::DaysFinderTool)
    optimize!(dft.m)
    stat = termination_status(dft.m)

    # Write results to dft if optimal
    if stat in [MOI.OPTIMAL, MOI.TIME_LIMIT] && has_values(dft.m)
        # Save selection variable as is
        dft.u = round.(getVariableValue(dft.misc[:u]))

        # Save representative periods and their mapping
        dft.rep_periods = [
            idx for (idx,val) in enumerate(dft.u) if val > 0
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

function getTimeSeriesErrorMatrix(dft::DaysFinderTool)
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
                sum(errorfunc.(
                    + dft.time_series[c].matrix_full_norm[p,t]
                    - dft.time_series[c].matrix_full_norm[pp,t]
                    for t in dft.timesteps
                ))
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

function getDurationCurveErrorMatrix(dft::DaysFinderTool)
    # TODO: This is a symmetric matrix, exploit that.
    squareBool = try_get_val(
        dft.config, "use_squared_error_matrix", false
    )
    errorfunc = squareBool ? x -> x^2 : x -> abs(x)
    DCEMD = Dict(
        c => [
            sum(
                errorfunc.(
                    .+ sort(dft.time_series[c].matrix_full_norm[p,:])
                    .- sort(dft.time_series[c].matrix_full_norm[pp,:])
                )
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

function getRandomHotStartForOrderingVariable(dft::DaysFinderTool)
    v_start_rand = spzeros(length.([dft.periods,dft.periods])...)
    idx = rand(dft.periods, dft.N_representative_periods - 1)
    v_start_rand[idx,idx] .= 1
    return v_start_rand
end

function fix_periods!(dft::DaysFinderTool)
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
