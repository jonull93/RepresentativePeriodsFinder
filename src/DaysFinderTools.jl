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

function runDaysFinderToolDefault(dft::DaysFinderTool, optimizer_factory)
    println("="^80)
    println("="^80)

    dft.AREA_ERROR_TOLERANCE = Dict()
    for c in dft.curves
        dft.AREA_ERROR_TOLERANCE[c] = dft.config["area_error_tolerance"]
    end

    u_val = Dict()
    w_val = Dict()
    for p in dft.periods
        u_val[p] = 0.
        w_val[p] = 0.
    end

    ###########################################################################
    # Model
    ###########################################################################
    m = Model(optimizer_factory)
    dft.m = m # save to DaysFinderTool

    ###########################################################################
    # Variables
    ###########################################################################
    @variable(m, u[p in dft.periods], Bin)
    @variable(m, w[p in dft.periods] >= 0)

    # Hot starts
    hot_start_values = try_get_val(dft.config, "hot_start_values", nothing)
    if isnothing(hot_start_values) == false
        u_val = hot_start_values[:u]
        w_val = hot_start_values[:w]
        for p in dft.periods
            set_start_value(u[p], u_val[p])
            set_start_value(w[p], w_val[p])
        end
    end

    dft.misc[:u] = u
    dft.misc[:w] = w

    # Mandatory Days
    m_d = []
    for ts in values(dft.time_series)
        mandatory_periods = get_mandatory_periods(ts, dft)
        append!(m_d, mandatory_periods)
        for p in mandatory_periods
            JuMP.fix(u[p],1)
            JuMP.set_lower_bound(w[p], 0.1)
        end
    end
    m_d = unique(m_d)

    # # Starting Values
    # p_start, w_start = define_random_start_point(dft, m_d, optimizer_factory)
    # for p in p_start
    #     @debug("Set starting value for $p")
    #     @debug("Set weight to $(w_start[p])")
    #     if ~is_fixed(u[p])
    #         set_start_value(u[p], 1)
    #     end
    #     set_start_value(w[p], w_start[p])
    # end

    ###########################################################################
    # Constraints
    ###########################################################################

    # Defining error as absolute value
    include_DC_error = try_get_val(dft.config, "include_DC_error", true)
    duration_curve_error = EmptyContainer(Float64)
    if include_DC_error == true
        @variable(m, duration_curve_error[c in dft.curves, b in dft.bins] >= 0)
        @constraint(m, error_eq1[c in dft.curves, b in dft.bins],
            duration_curve_error[c,b] >= + dft.L[c,b] - sum( w[p] / dft.N_total_periods * dft.A[c,p,b] for p in dft.periods)
        )
        @constraint(m, error_eq2[c in dft.curves, b in dft.bins],
            duration_curve_error[c,b] >= - dft.L[c,b] + sum( w[p] / dft.N_total_periods * dft.A[c,p,b] for p in dft.periods)
        )
    end

    # User defined number of representative periods
    @constraint(m, number_periods_eq,
        sum(u[p] for p in dft.periods) <= dft.N_representative_periods
    )

    # Restrict non-zero weights to selected periods
    if dft.config["equal_weights"] == true
        @constraint(m, single_weight_eq[p in dft.periods],
            w[p] == u[p] * dft.N_total_periods / dft.N_representative_periods
        )
    else
        @constraint(m, single_weight_eq[p in dft.periods],
            w[p] <= u[p] * dft.N_total_periods
        )
    end

    # Guarantee equivalent yearly duration
    @constraint(m, total_weight_eq,
        sum(w[p] for p in dft.periods) == dft.N_total_periods
    )

    # Minimum weight
    @constraint(m, minimum_weight[p in dft.periods],
        w[p] >= u[p] * 0.1
    )

    ##################################################################################
    # Optional constraints for speed up
    ##################################################################################
    # Defining area error
    enforce_area_error = try_get_val(dft.config, "enforce_area_error", true)
    if enforce_area_error == true
        @variable(m, area_error[c in dft.curves] >= 0)
        @constraint(m, area_error_eq1[c in dft.curves],
            area_error[c] >= dft.AREA_TOTAL[c] - sum(w[p] * dft.AREA_TOTAL_DAY[c,p] for p in dft.periods)
        )

        @constraint(m, area_error_eq2[c in dft.curves],
            area_error[c] <= dft.AREA_TOTAL[c] - sum(w[p] * dft.AREA_TOTAL_DAY[c,p] for p in dft.periods)
        )

        # Constraining area error
        @constraint(m, area_error_limit_eq[c in dft.curves],
            area_error[c] <= dft.AREA_ERROR_TOLERANCE[c]  * dft.AREA_TOTAL[c]
        )

        @constraint(m, helper1,
            sum(w[p] for p in dft.periods) >= sum(u[p] for p in dft.periods)
        )
    end

    ###########################################################################
    # Ordering of representative days constraints and variables
    ###########################################################################

    order_days = try_get_val(
        dft.config, "order_days", "none"
    )
    if order_days in ("binary","continuous")
        println("Timeseries error constraints...")
        # Define variable which maps days in the year to representative days
        # Can be continuous or binary
        if order_days == "continuous"
            @variable(m, 0 <= v[pp=dft.periods,p=dft.periods] <= 1)
        elseif order_days == "binary"
            @variable(m, v[pp=dft.periods,p=dft.periods], Bin)
        end
        dft.misc[:v] = v

        # Diagonal of v is equal to u
        @constraint(m, diag_chronology[p=dft.periods],
            v[p,p] == u[p]
        )

        # A day has to be assigned to at least one or multiple representative periods
        @constraint(m, [p=dft.periods],
            sum(v[pp,p] for pp=dft.periods) == 1
        )

        # The sum of the period has to be equal to the weight assigned to it
        @constraint(m, [p=dft.periods],
            sum(v[p,pp] for pp in dft.periods) == w[p]
        )

        # Helper constraint inspired by facility location problem
        @constraint(m, [p in dft.periods, pp in dft.periods],
            v[p,pp] <= u[p]
        )

        # # Define the reproduced timeseries
        # reproduced_timeseries = @expression(m,
        #     [c=dft.curves,p=dft.periods,t=dft.timesteps],
        #     sum(v[pp,p]*dft.time_series[c].matrix_full[pp,t] for pp=dft.periods)
        # )
        # dft.misc[:reproduced_timeseries] = reproduced_timeseries # save
        #
        # # Define the error for the timeseries
        # println("-"^80)
        # @variable(m, timeseries_error[c in dft.curves, p in dft.periods, t in dft.timesteps] >= 0)
        #
        # @constraint(m, [c in dft.curves, p in dft.periods, t in dft.timesteps],
        #     timeseries_error[c,p,t] >= reproduced_timeseries[c,p,t]
        #         - dft.time_series[c].matrix_full[p,t]
        # )
        #
        # @constraint(m, [c in dft.curves, p in dft.periods, t in dft.timesteps],
        #     timeseries_error[c,p,t] >= dft.time_series[c].matrix_full[p,t]
        #     - reproduced_timeseries[c,p,t]
        # )
    else
        v = EmptyContainer(Float64)
        timeseries_error = EmptyContainer(Float64)
    end

    ##################################################################################
    # Objective
    ##################################################################################

    TSEMD = getTimeSeriesErrorMatrix(dft)

    dc_weight = try_get_val(dft.config, "duration_curve_error_weight", 1.0)
    ts_weight = try_get_val(dft.config, "time_series_error_weight", 1.0)
    ramp_weight = try_get_val(dft.config, "ramping_error_weight", 1.0)

    @objective(m, Min,
        sum(
            dft.WEIGHT_DC[c] *(
                + dc_weight*sum(duration_curve_error[c,b] for b in dft.bins)
                + ts_weight*sum(
                    v[pp,p]*TSEMD[c][pp,p] for pp in dft.periods, p in dft.periods
                )
            ) for c in dft.curves
        )
    )

    optimize!(m)
    stat = termination_status(m)

    # Show results
    @show stat
    @show objective_value(m)
    @show objective_bound(m)
    @show has_values(m)

    if stat in [MOI.OPTIMAL, MOI.TIME_LIMIT] && has_values(m)
        u_val = value.(u).data
        w_val = value.(w).data
        dft.u = Dict(p => round(abs(u_val[p])) for p in dft.periods)
        dft.w = Dict(p => round(abs(w_val[p]), digits=4) for p in dft.periods)

        if order_days in ["binary", "continuous"]
            v_val = value.(v)
            dft.v = Dict(
                (i,j) => v_val[i,j] for i in v_val.axes[1], j in v_val.axes[2]
            )
        end
    else
        @show stat
    end
    return stat
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

    # If v[i,j] = 1, then day if is represented by day j
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
                    v[pp,p]*TSEMD[c][pp,p] for pp in dft.periods, p in dft.periods
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
    RP = [idx for (idx,val) in enumerate(u) if val[2] > 0] # selected periods j
    P = dft.periods # all periods i

    # Make ordering variable and weighting expression
    @variable(m, v[i=P,j=RP] >= 0)
    @expression(m, w[j=RP], sum(v[i,j] for i=P))
    dft.misc[:v] = v
    dft.misc[:w] = w

    # Enforce that representative period is selected for it's own period
    @constraint(m, [j=RP], v[j,j] == 1)

    # Ensure that every day is assigned a linear combination of rep periods
    @constraint(m, [i=P], sum(v[i,j] for j=RP) == 1)

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

##################################################################################
# Defining Starting Values
##################################################################################
function define_random_start_point(dft::DaysFinderTool, p_start::Array, optimizer_factory)
    w_start = Dict{String,Float64}()
    while true
        while length(p_start) < dft.N_representative_periods
            push!(p_start, rand(dft.periods))
            p_start = unique(p_start)
        end
        w_start = Dict()

        ##################################################################################
        # Optimization
        ##################################################################################
        m = Model(optimizer_factory)
        @variable(m, w[p in p_start] >= 0.1)
        @variable(m, area_error[c in dft.curves] >= 0)

        @objective(m, Min, sum(area_error[c] for c in dft.curves))

        # Guarantee equivalent yearly duration
        @constraint(m, total_weight_eq,
            sum(w[p] for p in p_start) == dft.N_total_periods
        )

        @constraint(m, area_error_eq1[c in dft.curves],
            area_error[c] >= + dft.AREA_TOTAL[c] - sum(w[p] * dft.AREA_TOTAL_DAY[c,p] for p in p_start)
        )

        @constraint(m, area_error_eq2[c in dft.curves],
            area_error[c] <= + dft.AREA_TOTAL[c] - sum(w[p] * dft.AREA_TOTAL_DAY[c,p] for p in p_start)
        )

        # Constraining area error
        @constraint(m, area_error_limit_eq[c in dft.curves],
            area_error[c] <= dft.AREA_ERROR_TOLERANCE[c]  * dft.AREA_TOTAL[c]
        )

        # print(m)
        optimize!(m)
        stat = termination_status(m)
        for k in keys(w)
            w_start[k[1]] = round(abs(value(w[k[1]])), digits=4)
        end

        @debug("Status finding starting values: $stat")

        if stat in [MOI.OPTIMAL]
            break
        end
        for c in dft.curves
            dft.AREA_ERROR_TOLERANCE[c] += 0.01
        end
    end
    @debug("$w_start")
    return p_start, w_start
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
