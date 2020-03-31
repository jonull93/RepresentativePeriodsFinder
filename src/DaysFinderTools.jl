##################################################################################
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
    periods::UnitRange                # set of potential representative periods
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
    u::Dict
    w::Dict
    v::Dict

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
    # Sets
    ##################################################################################
    self.bins =     ["b"*string(b, pad=pad) for b in range(1,stop=self.config["number_bins"])]

    lenT = try_get_val(self.config, "timesteps_per_period", 24)
    self.timesteps = 1:lenT

    lenP = try_get_val(self.config, "number_days_total", 365)
    self.periods = 1:lenP


    ##################################################################################
    # Parameters
    ##################################################################################
    self.WEIGHT_DC                  = Dict()
    self.A                          = Dict()
    self.L                          = Dict()

    self.AREA_TOTAL                 = Dict()
    self.AREA_TOTAL_DAY             = Dict()

    self.N_representative_periods   = self.config["number_days"]
    self.N_total_periods            = self.config["number_days_total"]

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
    @variable(m, w[j in dft.periods] >= 0)
    @variable(m, v[i in dft.periods, j in dft.periods] >= 0)
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
        sum(u[j] for j in dft.periods) <= dft.N_representative_periods
    )

    # Restrict non-zero weights to selected periods
    equal_weights = try_get_val(dft.config, "equal_weights", true)
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

    # Constraints for v - no effect on results, just for things to work smoothly
    # If day j is selected as a representative day then u[j] is 1
    @constraint(m, [i in dft.periods, j in dft.periods],
        v[i,j] <= u[j]
    )
    @constraint(m, [j in dft.periods],
        sum(v[i,j] for i in dft.periods) <= w[j]
    )

    # Objective
    dc_weight = try_get_val(dft.config, "duration_curve_error_weight", 1.0)
    obj = @expression(m,
        sum(
            dft.WEIGHT_DC[c] *(
                + dc_weight*sum(duration_curve_error[c,b] for b in dft.bins)
            ) for c in dft.curves
        )
    )
    @objective(m, Min, obj)
    return m
end

function makeDaysFinderToolModel(dft::DaysFinderTool, optimizer_factory)
    # Model
    m = Model(optimizer_factory)
    dft.m = m

    # Create variables
    @variable(m, u[j in dft.periods], Bin)
    @variable(m, v[i in dft.periods, j in dft.periods], Bin)
    @variable(m, w[j in dft.periods])
    @constraint(m, [j in dft.periods],
        sum(v[i,j] for i in dft.periods) == w[j]
    )
    dft.misc[:u] = u
    dft.misc[:v] = v
    dft.misc[:w] = w

    # Hot start v
    v_start = try_get_val(
        dft.config, "hot_start_values",
        Dict((i,j) => i < dft.N_total_periods ? 1.0 : 0.0 for i in dft.periods, j in dft.periods)
    )
    for i in dft.periods, j in dft.periods
        set_start_value(v[i,j], v_start[i,j])
    end

    # If v[i,j] = 1, then day i is represented by day j
    # The sum over j (columns) is the weighting of a particular day
    # The sum over i (rows) must equal 1 - a day is represented by only one other day

    # Each day i has to be represented by another day j
    @constraint(m, [i in dft.periods],
        sum(v[i,j] for j in dft.periods) == 1
    )

    # If day j is selected as a representative day then u[j] is 1
    @constraint(m, [i in dft.periods, j in dft.periods],
        v[i,j] <= u[j]
    )

    # Choose only N_representative periods
    @constraint(m,
        sum(u[j] for j in dft.periods) <= dft.N_representative_periods
    )

    # Sum of weightings of representative periods must equal N_total_periods
    @constraint(m,
        sum(w[j] for j in dft.periods) == dft.N_total_periods
    )

    # Define DC error
    @variable(m, duration_curve_error[c in dft.curves, b in dft.bins] >= 0)
    @constraint(m, error_eq1[c in dft.curves, b in dft.bins],
        duration_curve_error[c,b] >= + dft.L[c,b] - sum(w[j] / dft.N_total_periods * dft.A[c,j,b] for j in dft.periods)
    )
    @constraint(m, error_eq2[c in dft.curves, b in dft.bins],
        duration_curve_error[c,b] >= - dft.L[c,b] + sum(w[j] / dft.N_total_periods * dft.A[c,j,b] for j in dft.periods)
    )

    # Other constraints
    # TODO Impose particular day as a solution

    # Define objective
    TSEMD = getTimeSeriesErrorMatrix(dft)
    dc_weight = try_get_val(dft.config, "duration_curve_error_weight", 1.0)
    ts_weight = try_get_val(dft.config, "time_series_error_weight", 1.0)
    obj = @expression(m,
        sum(
            dft.WEIGHT_DC[c] *(
                + dc_weight*sum(duration_curve_error[c,b] for b in dft.bins)
                + ts_weight*sum(
                    v[pp,j]*TSEMD[c][pp,p] for pp in dft.periods, p in dft.periods
                )
            ) for c in dft.curves
        )
    )
    @objective(m, Min, obj)

    # Constrain lower bound of objective
    obj_lower_bound = try_get_val(
        dft.config, "objective_lower_bound", 0.0
    )
    @constraint(m, obj >= obj_lower_bound)

    return m
end

function makeReOrderingDaysFinderTool(dft::DaysFinderTool, optimizer_factory)
    println("-" ^ 80)
    println("-" ^ 80)
    println("Re ordering days")
    println("-" ^ 80)
    println("-" ^ 80)

    # Make model
    m = Model(optimizer_factory)
    dft.m = m

    # Retrieve selection variable
    u = dft.u

    # Get sets
    RP = [val[1] for (idx,val) in enumerate(u) if val[2] > 0] # selected periods j
    P = dft.periods # all periods i

    # Make ordering variable and weighting expression
    @variable(m, v[i=P,j=RP] >= 0)
    @expression(m, w[j=RP], sum(v[i,j] for i=P))
    dft.misc[:v] = v
    dft.misc[:w] = w

    # Enforce same weightings if applicable
    enforce_same_weightings = try_get_val(
        dft.config, "enforce_same_weightings", false
    )
    if enforce_same_weightings == true
        @constraint(m, [j=RP],
            w[j] == dft.w[j]
        )
    end

    # Enforce that representative period is selected for it's own period
    @constraint(m, [j=RP], v[j,j] == 1)

    # Ensure that every day is assigned a linear combination of rep periods
    linear_comb = try_get_val(
        dft.config["solver"], "LinearCombination", "sum_to_one"
    )
    if linear_comb == "sum_to_1"
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
        sum(v[i,j]*dft.time_series[c].matrix_full[j,t] for j=RP)
    )

    # Define error
    @variable(m, ts_error[c=dft.curves, i=P, t=dft.timesteps] >= 0)
    @constraint(m,
    	[c=dft.curves, i=P, t=dft.timesteps],
    	ts_error[c,i,t] >=
            + reconstructed_time_series[c,i,t]
    		- dft.time_series[c].matrix_full[i,t]
    )
    @constraint(m,
    	[c=dft.curves, i=P,t=dft.timesteps],
    	ts_error[c,i,t] >= -(
            + reconstructed_time_series[c,i,t]
    		- dft.time_series[c].matrix_full[i,t]
        )
    )

    # Define objective
    TSEMD = getTimeSeriesErrorMatrix(dft)
    dc_weight = try_get_val(dft.config, "duration_curve_error_weight", 1.0)
    ts_weight = try_get_val(dft.config, "time_series_error_weight", 1.0)

    @objective(m, Min,
        sum(
            dft.WEIGHT_DC[c] *(
                + dc_weight*sum(duration_curve_error[c,b] for b in dft.bins)
                + ts_weight*sum(ts_error[c,i,t] for i in P, t in dft.timesteps)
            ) for c in dft.curves
        )
    )

    return m
end

function optimizeDaysFinderTool(dft::DaysFinderTool)
    optimize!(dft.m)
    stat = termination_status(dft.m)

    # Write results to dft if optimal
    if stat in [MOI.OPTIMAL, MOI.TIME_LIMIT] && has_values(dft.m)
        u_val = value.(dft.misc[:u])
        v_val = value.(dft.misc[:v])
        w_val = value.(dft.misc[:w])

        dft.u = Dict(
            j => j in u_val.axes[1] ? u_val[j] : 0 for j in dft.periods
        )
        dft.w = Dict(
            j => j in w_val.axes[1] ? w_val[j] : 0 for j in dft.periods
        )
        dft.v = Dict(
            (i,j) => j in v_val.axes[2] ? v_val[i,j] : 0 for i in dft.periods, j in dft.periods
        )
    else
        @show stat
    end

    # Return solution status
    return stat
end

function getTimeSeriesErrorMatrix(dft::DaysFinderTool)
    TSEM = Array{Float64,2}(undef, length(dft.periods), length(dft.periods))
    TSEMD = Dict(
        c => abs.([sum(
            + dft.time_series[c].matrix_full[p,t]
            - dft.time_series[c].matrix_full[pp,t] for t in dft.timesteps
        ) for  p in dft.periods, pp in dft.periods]) for c in dft.curves
    )
    return TSEMD
end
