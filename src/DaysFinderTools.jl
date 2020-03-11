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

##################################################################################
#  Add TimeSeries to Model
##################################################################################
function addTimeSeries!(self::DaysFinderTool, ts::TimeSeries)
    self.time_series[ts.name] = ts
    push!(self.curves, ts.name)

    weight!(self.WEIGHT_DC, ts)
    area_total!(self.AREA_TOTAL, ts)
    area_total_days!(self.AREA_TOTAL_DAY, self.periods, ts)

    cum_bin_end!(self.A, self.periods, self.bins, ts)
    cum_bin_total!(self.L, self.bins, ts)
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

##################################################################################
# Write out results to  csv and yaml files
##################################################################################
function writeOutResults(dft::DaysFinderTool)
    result_dir = normpath(joinpath(dft.config["basedir"], dft.config["result_dir"]))
    if !isdir(result_dir)
        mkpath(result_dir)
    end

    df_dv = DataFrame(
                periods     = sort([k for k in dft.periods]),
                weights     = [dft.w[k] for k in sort([k for k in keys(dft.w)])],
                used_days   = [dft.u[k] for k in sort([k for k in keys(dft.u)])])
    CSV.write(joinpath(result_dir, "decision_variables.csv"), df_dv, delim=';')

    df_dv_s = deepcopy(df_dv[df_dv[!,:weights] .> 0, :])
    CSV.write(joinpath(result_dir, "decision_variables_short.csv"), df_dv_s, delim=';')

    ###########################################################################
    # Resulting time series
    ###########################################################################
    df_ts = DataFrame()
    period_idx = df_dv[!,:used_days] .== 1
    @debug("Period index of selected days: $period_idx")

    for ts in values(dft.time_series)
        df_ts[!,Symbol(ts.name)] =  ts.matrix_full[period_idx,:]'[:]
    end
    CSV.write(joinpath(result_dir, "resulting_profiles.csv"), df_ts, delim=';')

    ###########################################################################
    # Save the ordering variable (if necessary)
    ###########################################################################
    order_days = try_get_val(
        dft.config, "order_days", false
    )
    if order_days
        IPW = [dft.v[pp,p] for p in dft.periods, pp in dft.periods]
        df = DataFrame(IPW)
        CSV.write(joinpath(result_dir, "ordering_variable.csv"), df, delim=';')
    end

    ###########################################################################
    # Copy of config-file
    ###########################################################################
    stringdata = JSON.json(dft.config)
    open(joinpath(result_dir, "config_file.json"), "w") do f
        write(f, stringdata)
    end
    return dft
end

###############################################################################
# Run Optimization
###############################################################################
# import RepresentativeDaysFinders: runDaysFinderToolDefault
# using RepresentativeDaysFinders: DaysFinderTool
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
    @variable(m, area_error[c in dft.curves] >= 0)

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

    duration_curve_error = try_get_val(
        dft.config, "duration_curve_error", "absolute"
    )
    if duration_curve_error == "absolute"
        # Defining error as absolute value
        @variable(m, duration_curve_error[c in dft.curves, b in dft.bins] >= 0)
        @constraint(m, error_eq1[c in dft.curves, b in dft.bins],
            duration_curve_error[c,b] >= + dft.L[c,b] - sum( w[p] / dft.N_total_periods * dft.A[c,p,b] for p in dft.periods)
        )
        @constraint(m, error_eq2[c in dft.curves, b in dft.bins],
            duration_curve_error[c,b] >= - dft.L[c,b] + sum( w[p] / dft.N_total_periods * dft.A[c,p,b] for p in dft.periods)
        )
    elseif duration_curve_error == "square"
        # Define error as 2 norm error
        duration_curve_error = Dict()
        for c in dft.curves, b in dft.bins
            duration_curve_error[c,b] = @expression(m,
                (
                    dft.L[c,b] - sum(w[p] / dft.N_total_periods * dft.A[c,p,b] for p in dft.periods)
                )^2
            )
        end
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

    ###########################################################################
    # Ordering of representative days constraints and variables
    ###########################################################################

    order_days = try_get_val(
        dft.config, "order_days", false
    )
    if order_days == true

        println("Timeseries error constraints...")

        # Define variable which maps days in the year to representative days
        # Can be continuous or binary
        # TODO: allow for the option
        if true
            @variable(m, 0 <= v[pp=dft.periods,p=dft.periods] <= 1)
        elseif false
            @variable(m, v[pp=dft.periods,p=dft.periods], Bin)
        else
            error("Please specify whether ordering variable should be binary or continuous")
        end

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

        # Define the reproduced timeseries
        reproduced_timeseries = @expression(m,
            [c=dft.curves,p=dft.periods,t=dft.timesteps],
            sum(v[pp,p]*dft.time_series[c].matrix_full[pp,t] for pp=dft.periods)
        )
        dft.misc[:reproduced_timeseries] = reproduced_timeseries # save

        # Define the error for the timeseries
        println("-"^80)
        @variable(m, timeseries_error[c in dft.curves, p in dft.periods, t in dft.timesteps] >= 0)

        @constraint(m, [c in dft.curves, p in dft.periods, t in dft.timesteps],
            timeseries_error[c,p,t] >= reproduced_timeseries[c,p,t]
                - dft.time_series[c].matrix_full[p,t]
        )

        @constraint(m, [c in dft.curves, p in dft.periods, t in dft.timesteps],
            timeseries_error[c,p,t] >= dft.time_series[c].matrix_full[p,t]
            - reproduced_timeseries[c,p,t]
        )

        # c = dft.curves[1]; p = dft.periods[1]; t = dft.timesteps[1]
        # @show reproduced_timeseries[c,p,t]
        # @show dft.time_series[c].matrix_full[p,t]
        # timeseries_error = Dict(
        #     (c,p,t) => @expression(m,
        #         (dft.time_series[c].matrix_full[p,t]
        #         - reproduced_timeseries[c,p,t])^2
        #     ) for c in dft.curves, p in dft.periods, t in dft.timesteps
		# )

        # timeseries_error = @expression(m,
        #     [c=dft.curves,p=dft.periods,t=dft.timesteps],
        #     dft.time_series[c].matrix_full[p,t]
        #     - reproduced_timeseries[c,p,t]
        # )

        # timeseries_error = Dict(
        #     (c,p,t) => 0.0 for c in dft.curves, p in dft.periods, t=dft.timesteps
        # )
    else
        timeseries_error = Dict(
            (c,p,t) => 0.0 for c in dft.curves, p in dft.periods, t=dft.timesteps
        )
    end

    # Save some shit to the misc dictionary
    dft.misc[:u] = u
    order_days ? dft.misc[:v] = v : dft.misc[:v] = Dict()

    ##################################################################################
    # Objective
    ##################################################################################

    dc_weight = try_get_val(dft.config, "duration_curve_error_weight", 1.0)
    ts_weight = try_get_val(dft.config, "time_series_error_weight", 1.0)
    ramp_weight = try_get_val(dft.config, "ramping_error_weight", 1.0)

    @objective(m, Min,
        sum(
            dft.WEIGHT_DC[c] *(
                + dc_weight*sum(duration_curve_error[c,b] for b in dft.bins)
                + ts_weight*sum(timeseries_error[c,p,t] for p in dft.periods, t=dft.timesteps)
            ) for c in dft.curves
        )
    )

    # @objective(m, Min,
    #     sum(
    #         dft.WEIGHT_DC[c] *(
    #             + dc_weight*sum(dft.L[c,b] - sum(w[p] / dft.N_total_periods * dft.A[c,p,b] for p in dft.periods) for b in dft.bins)
    #             + ts_weight*sum( for p in dft.periods, b in dft.bins)
    #         ) for c in dft.curves
    #     )
    # )

    # @objective(m, Min, sum(dft.WEIGHT_DC[c]*
    #     (dft.time_series[c].matrix_full[p,t] - reproduced_timeseries[c,p,t])^2
    #     for c in dft.curves, p in dft.periods, t in dft.timesteps
    # ))

    optimize!(m)
    stat = termination_status(m)

    @show stat
    @show objective_value(m)
    @show objective_bound(m)
    @show has_values(m)

    if stat in [MOI.OPTIMAL, MOI.TIME_LIMIT] && has_values(m)
        u_val = value.(u)
        w_val = value.(w)

        dft.u = Dict()
        for k in keys(u)
            dft.u[k[1]] = round(Int,abs(value(u[k[1]])))
        end
        dft.w = Dict()
        for k in keys(w)
            dft.w[k[1]] = round(abs(value(w[k[1]])), digits=4)
        end
        if order_days == true
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
