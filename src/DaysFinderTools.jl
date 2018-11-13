##################################################################################
# Author:   Hanspeter Höschle
# Date:     15/06/2017
##################################################################################
type DaysFinderTool

    ##################################################################################
    # Configuration
    ##################################################################################
    config_file::AbstractString
    config::Dict{Any,Any}

    ##################################################################################
    # Sets
    ##################################################################################
    bins::Array{String,1}                   # set of bins
    curves::Array{String,1}                 # set of duration curves
    periods::Array{String,1}                # set of potential representative periods

    ##################################################################################
    # Parameters
    ##################################################################################
    N_representative_periods::Int           # Number of representative periods to select
    N_total_periods::Int                    # Total number of repetitions required to scale up the duration of a single representative period to one year

    WEIGHT_DC::Dict{String,Float64}                         # Relative importance of approximating duration curve c
    A::Dict{Tuple{String,String,String},Float64}            # Share of the time of day d during which the lowest value of the range corresponding to bin b of duration curve c is exceeded
    L::Dict{Tuple{String,String},Float64}                   # Share of the time during which the values of a time series with corresponding duration curve c exceed the lowest value of the range corresponding to bin b
    AREA_TOTAL::Dict{String,Float64}                        # Total area under the DC of time series c
    AREA_TOTAL_DAY::Dict{Tuple{String,String},Float64}      # Total area under time series with DC c per day
    AREA_ERROR_TOLERANCE::Dict{String,Float64}              # tolerance parameter for area constraint

    ##################################################################################
    # Time Series
    ##################################################################################
    time_series::Dict{String,TimeSeries}

    ##################################################################################
    # Results
    ##################################################################################
    u::Dict{String,Int}
    w::Dict{String,Float64}

    function DaysFinderTool(config_file::String)
        self = new()
        self.config_file = config_file
        self.config = YAML.load(open(config_file))

        ##################################################################################
        # Sets
        ##################################################################################
        self.bins =     [@sprintf "b%03d" b for b in range(1,self.config["number_bins"])]
        self.periods =  [@sprintf "p%03d" p for p in range(1,self.config["number_days_total"])]

        ##################################################################################
        # Parameters
        ##################################################################################
        self.WEIGHT_DC                  = Dict{String,Float64}()
        self.A                          = Dict{Tuple{String,String,String},Float64}()
        self.L                          = Dict{Tuple{String,String},Float64}()

        self.AREA_TOTAL                 = Dict{String,Float64}()
        self.AREA_TOTAL_DAY             = Dict{Tuple{String,String},Float64}()

        self.N_representative_periods   = self.config["number_days"]
        self.N_total_periods            = self.config["number_days_total"]

        ##################################################################################
        # Initialize the TimeSeries
        ##################################################################################
        self.time_series = Dict{String,TimeSeries}()
        self.curves = Array{String,1}()
        for ts_config in self.config["time_series"]
            @printf "%s\n" repeat("-", 100)

            ts = TimeSeries(self, ts_config)
            addTimeSeries!(self, ts)

            @printf "%-20s added\n" ts.name

            ##################################################################################
            # Dynamics profiles
            ##################################################################################
            if self.config["dynamics_method"] == "all_1step"
                @printf "%s\n" repeat("-", 100)

                ts_diff = TimeSeries(self, ts, :all_1step)
                addTimeSeries!(self, ts_diff)

                @printf "%-20s added\n" ts_diff.name
            end
        end

        ##################################################################################
        # Correlation
        ##################################################################################
        if self.config["correlation_method"] == "all"
            ts_basic = [ts.name for ts in values(self.time_series) if ts.time_series_type ==:basic]
            for combi in combinations(ts_basic,2)
                @printf "%s\n" repeat("-", 100)
                ts1 = self.time_series[combi[1]]
                ts2 = self.time_series[combi[2]]
                ts = TimeSeries(self, ts1, ts2)
                addTimeSeries!(self, ts)

                @printf "%-20s added\n" ts.name
            end
        end
        return self
    end
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
# Run Optimization
##################################################################################
function runDaysFinderToolIterativeBounderaries(dft::DaysFinderTool)
    @printf "%s\n" repeat("=", 100)
    @printf "%s\n" repeat("=", 100)

    u_val = Dict{String,Int}()
    w_val = Dict{String,Int}()
    for p in dft.periods
        u_val[p] = 0.
        w_val[p] = 0.
    end
    Big_M_wp_step = (dft.N_total_periods - (dft.N_total_periods / dft.N_representative_periods)) / 50.
    Big_M_wp = dft.N_total_periods / dft.N_representative_periods + Big_M_wp_step

    log("debug", "Step increase Big-M: $Big_M_wp_step")
    TOL_MIN = Dict{String,Float64}()
    for c in dft.curves
        if contains(c, "diff")
            TOL_MIN[c] = 1e-4
        elseif contains(c, "corr")
            TOL_MIN[c] = 1e-4
        else
            TOL_MIN[c] =  1e-8
        end
    end
    TIME_LIMIT = length(dft.curves) * 2.5
    log("debug", "Time limit for iterative runs: $TIME_LIMIT")
    LONG_RUN = false

    stat = 0

    while true

        log("debug", "Big-M for weight constraint: $Big_M_wp")
        log("debug", "Min-Tolerance for user cuts: $TOL_MIN")

        if ~LONG_RUN
            dft.AREA_ERROR_TOLERANCE = Dict{String,Float64}()
            for c in dft.curves
                dft.AREA_ERROR_TOLERANCE[c] = dft.config["area_error_tolerance"]
            end
        end

        ##################################################################################
        # Model
        ##################################################################################
        # if dft.config["solver"]["mip"] == "CPLEX"
            # m = Model(
            #     solver=CplexSolver(
            #                 CPX_PARAM_SCRIND=1,
            #                 CPX_PARAM_TILIM=TIME_LIMIT,
            #                 CPX_PARAM_THREADS=dft.config["solver"]["Threads"]))

        if dft.config["solver"]["mip"] == "GLPK"
            m = Model(solver=GLPKSolverMIP(presolve=true, msg_lev=GLPK.MSG_ON))
            # getSolverModel(m).params.tm_lim=15


        elseif dft.config["solver"]["mip"] == "CBC"
            m = Model(
                    solver=CbcSolver(
                                seconds=dft.config["solver"]["ResLim"],
                                logLevel=1,
                                threads=32))

        elseif dft.config["solver"]["mip"] == "Gurobi"
            m = Model(
                    solver=GurobiSolver(
                        Presolve=0,
                        TimeLimit=TIME_LIMIT,
                        Threads=dft.config["solver"]["Threads"]
                        )
                    )

        else
            error()
            return
        end
        ##################################################################################
        # Variables
        ##################################################################################
        @variable(m, u[p in dft.periods], Bin, start=u_val[p])
        @variable(m, w[p in dft.periods] >= 0, start=w_val[p])

        @variable(m, area_error[c in dft.curves] >= 0)
        @variable(m, error[c in dft.curves, b in dft.bins] >= 0)

        # Mandatory Days
        m_d = String[]
        for ts in values(dft.time_series)
            mandatory_periods = get_mandatory_periods(ts, dft)
            append!(m_d, mandatory_periods)
            for p in mandatory_periods
                JuMP.fix(u[p],1)
                setlowerbound(w[p], 0.1)
            end
        end
        m_d = unique(m_d)

        # Starting Values
        # p_start, w_start = define_random_start_point(dft, m_d)
        # for p in p_start
        #     log("debug", "Set starting value for $p")
        #     log("debug", "Set weight to $(w_start[p])")
        #     if getcategory(u[p]) != :Fixed
        #         setvalue(u[p], 1)
        #     end
        #     setvalue(w[p], w_start[p])
        # end

        ##################################################################################
        # Objective
        ##################################################################################
        @objective(m, Min,
            + sum(dft.WEIGHT_DC[c] * sum(error[c,b] for b in dft.bins) for c in dft.curves)
            + 1 * sum(area_error[c] for c in dft.curves) # to bring area_error to lower bound
        )

        ##################################################################################
        # Constraints
        ##################################################################################
        # Defining error as absolute value
        @constraint(m, error_eq1[c in dft.curves, b in dft.bins],
            error[c,b] >= + dft.L[c,b] - sum( w[p] / dft.N_total_periods * dft.A[c,p,b] for p in dft.periods)
        )
        @constraint(m, error_eq2[c in dft.curves, b in dft.bins],
            error[c,b] >= - dft.L[c,b] + sum( w[p] / dft.N_total_periods * dft.A[c,p,b] for p in dft.periods)
        )

        # User defined number of representative periods
        @constraint(m, number_periods_eq,
            sum(u[p] for p in dft.periods) <= dft.N_representative_periods
        )

        # Restrict non-zero weights to selected periods
        @constraint(m, single_weight_eq[p in dft.periods],
            w[p] <= u[p] * Big_M_wp
        )

        # Guarantee equivalent yearly duration
        @constraint(m, total_weight_eq,
            sum(w[p] for p in dft.periods) == dft.N_total_periods
        )

        # Minimum weight
        log("debug","minimum weight set to 0.1")
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
            area_error[c] >= - dft.AREA_TOTAL[c] + sum(w[p] * dft.AREA_TOTAL_DAY[c,p] for p in dft.periods)
        )

        # Constraining area error
        @constraint(m, area_error_limit_eq[c in dft.curves],
            area_error[c] <= dft.AREA_ERROR_TOLERANCE[c]  * dft.AREA_TOTAL[c]
        )

        @constraint(m, helper1,
            sum(w[p] for p in dft.periods) >= sum(u[p] for p in dft.periods)
        )
        ##################################################################################
        # Callback-function
        ##################################################################################

        function mycutgenerator(cb)
            area_error_val = getvalue(area_error)

            for c in dft.curves
                # @sprintf "%s" repeat("=", 80)
                # @show c
                # @show area_error_val[c]
                # @show dft.AREA_TOTAL[c]
                NEW_AREA_ERROR_TOLERANCE = max(abs(area_error_val[c] / dft.AREA_TOTAL[c]), TOL_MIN[c])

                # @show NEW_AREA_ERROR_TOLERANCE
                # @show dft.AREA_ERROR_TOLERANCE

                if NEW_AREA_ERROR_TOLERANCE < dft.AREA_ERROR_TOLERANCE[c]
                    log("info", "Add new area error cut for $c ")
                    log("debug", "$NEW_AREA_ERROR_TOLERANCE, $(dft.AREA_ERROR_TOLERANCE[c]) ")
                    @usercut(cb, area_error[c] <= NEW_AREA_ERROR_TOLERANCE * dft.AREA_TOTAL[c])
                    dft.AREA_ERROR_TOLERANCE[c] = NEW_AREA_ERROR_TOLERANCE
                end
            end
        end
        if ~LONG_RUN
            addcutcallback(m, mycutgenerator)
        end

    # print(m)

        stat = solve(m)
        u_val = getvalue(u)
        w_val = getvalue(w)

        area_error_val = getvalue(area_error)
        area_error_real_1 = Dict{String, Float64}()
        area_error_real_2 = Dict{String, Float64}()
        area_error_tol = Dict{String, Float64}()

        for c in dft.curves
            sum_area = sum([w_val[p] .* dft.AREA_TOTAL_DAY[c,p] for p in dft.periods])
            # @show sum_area
            area_error_real_1[c] = dft.AREA_TOTAL[c] - sum_area
            area_error_real_2[c] = - dft.AREA_TOTAL[c] + sum_area
            area_error_tol[c] = dft.AREA_ERROR_TOLERANCE[c]  * dft.AREA_TOTAL[c]
        end
        for c in dft.curves
            log("debug", "$c: $(area_error_real_1[c]) $(area_error_real_2[c]) <=  $(area_error_val[c]) <=  $(area_error_tol[c])")
        end

        test_Big_M = [w_val[p] < u_val[p] * Big_M_wp - 1e-6 for p in dft.periods if u_val[p] == 1]
        test_Area_Error = [area_error_val[c] < dft.AREA_ERROR_TOLERANCE[c]  * dft.AREA_TOTAL[c] for c in dft.curves]

        res = [(p,w_val[p], Big_M_wp - w_val[p]) for p in dft.periods if u_val[p] == 1]
        @show  test_Big_M all(test_Big_M)  (length(test_Big_M) - count(test_Big_M))
        @show  test_Area_Error all(test_Area_Error) (length(test_Area_Error) - count(test_Area_Error))
        @show Big_M_wp res

        if all(test_Big_M) && all(test_Area_Error)
            if TIME_LIMIT == dft.config["solver"]["ResLim"]
                break
            else
                TIME_LIMIT = dft.config["solver"]["ResLim"]
                LONG_RUN = true
            end
        elseif isempty(test_Big_M)
            Big_M_wp  += Big_M_wp_step
            LONG_RUN = false
        elseif all(test_Big_M) && ~all(test_Area_Error)
            for (i,c) in enumerate(dft.curves)
                if ~test_Area_Error[i]
                    TOL_MIN[c] *= 10
                end
            end
            TIME_LIMIT = length(dft.curves) * 2.5
            LONG_RUN = false
        else
            Big_M_wp  += Big_M_wp_step
            TIME_LIMIT = length(dft.curves) * 2.5
            LONG_RUN = false
        end
    end


    if stat in [:Optimal, :UserLimit]
        u = getvalue(u)
        dft.u = Dict{String,Int}()
        for k in keys(u)
            dft.u[k[1]] = round(Int,abs(u[k[1]]))
        end

        w = getvalue(w)
        dft.w = Dict{String,Int}()
        for k in keys(w)
            dft.w[k[1]] = round(abs(w[k[1]]),4)
        end

    else
        @show stat
    end
    return stat
end

##################################################################################
# Brute Force
##################################################################################
function runDaysFinderToolBruteForce(dft::DaysFinderTool)
    TIME_LIMIT = dft.config["solver"]["ResLim"]
    ITERATION_LIMIT = dft.config["solver"]["ITERATION_LIMIT"]
    iteration = 1

    ##################################################################################
    # Results
    ##################################################################################
    result_obj = 1e12
    result_u = 0
    result_w = 0

    u_val_zeros = Dict{String,Float64}()
    for d in dft.periods
        u_val_zeros[d] = 0
    end

    # Mandatory Periods
    mandatory_periods = String[]
    for ts in values(dft.time_series)
        m_p = get_mandatory_periods(ts, dft)
        append!(mandatory_periods, m_p)
    end
    mandatory_periods = unique(mandatory_periods)
    for p in mandatory_periods
        u_val_zeros[p] = 1
    end


    time = 0
    while true
        tic()
        log("info", "Iteration: $iteration")
        u_val = copy(u_val_zeros)
        while sum(values(u_val)) < dft.N_representative_periods
            u_val[rand(dft.periods,1)[1]] = 1
        end
        # @show sum(values(u_val))
        ##################################################################################
        # Model
        ##################################################################################
        # if dft.config["solver"]["mip"] == "CPLEX"
        #     m = Model(solver=CplexSolver(CPX_PARAM_SCRIND=0))
        if dft.config["solver"]["mip"] == "Gurobi"
            m = Model(solver=GurobiSolver(OutputFlag=0))
        end

        @variable(m, w[p in dft.periods] >= 0)
        @variable(m, u[p in dft.periods], Bin)
        @variable(m, error[c in dft.curves, b in dft.bins] >= 0)
        for p in mandatory_periods
            setlowerbound(w[p], 0.1)
        end

        ##################################################################################
        # Objective
        ##################################################################################
        @objective(m, Min, sum(dft.WEIGHT_DC[c] * sum(error[c,b] for b in dft.bins) for c in dft.curves))

        ##################################################################################
        # Constraints
        ##################################################################################
        # Defining error as absolute value
        @constraint(m, error_eq1[c in dft.curves, b in dft.bins],
            error[c,b] >= + dft.L[c,b] - sum( w[p] / dft.N_total_periods * dft.A[c,p,b] for p in dft.periods)
        )
        @constraint(m, error_eq2[c in dft.curves, b in dft.bins],
            error[c,b] >= - dft.L[c,b] + sum( w[p] / dft.N_total_periods * dft.A[c,p,b] for p in dft.periods)
        )

        # User defined number of representative periods
        @constraint(m, number_periods_eq,
            sum(u_val[p] for p in dft.periods) == dft.N_representative_periods
        )

        # Restrict non-zero weights to selected periods
        @constraint(m, single_weight_eq[p in dft.periods],
            w[p] <= u_val[p] * dft.N_total_periods
        )

        # Guarantee equivalent yearly duration
        @constraint(m, total_weight_eq,
            sum(w[p] for p in dft.periods) == dft.N_total_periods
        )

        @constraint(m, helper[p in dft.periods],
            u[p] == u_val[p])

        stat = solve(m)
        obj = getobjectivevalue(m)
        # @show stat
        if stat == :Optimal && obj < result_obj
            result_u = getvalue(u)
            result_w = getvalue(w)
            result_obj = obj
            log("info", "found improved objective: $result_obj")
        end
        time += toq()
        # @show time
        if (iteration >= ITERATION_LIMIT) | (time >= TIME_LIMIT)
            break
        end
        iteration += 1
    end

    log("info", "Best found objective: $result_obj after $iteration Iterations and $time seconds")

    dft.u = Dict{String,Int}()
    for k in keys(result_u)
        dft.u[k[1]] = round(Int,abs(result_u[k[1]]))
    end

    dft.w = Dict{String,Int}()
    for k in keys(result_w)
        dft.w[k[1]] = round(abs(result_w[k[1]]),4)
    end

    return :UserLimit
end


##################################################################################
# Defining Starting Values
##################################################################################
function define_random_start_point(dft::DaysFinderTool, p_start::Array{String})
    w_start = Dict{String,Float64}()
    while true
        while length(p_start) < dft.N_representative_periods
            push!(p_start, rand(dft.periods))
            p_start = unique(p_start)
        end
        w_start = Dict{String,Float64}()

        ##################################################################################
        # Optimization
        ##################################################################################
        m = Model(solver=GurobiSolver())
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
        stat = solve(m)
        w = getvalue(w)
        for k in keys(w)
            w_start[k[1]] = round(abs(w[k[1]]),4)
        end

        log("debug", "Status finding starting values: $stat")

        if stat in [:Optimal]
            break
        end
        for c in dft.curves
            dft.AREA_ERROR_TOLERANCE[c] += 0.01
        end
    end
    log("debug", "$p_start")
    log("debug", "$w_start")
    return p_start, w_start
end

##################################################################################
# Write out results to  csv and yaml files
##################################################################################
function writeOutResults(dft::DaysFinderTool)
    result_dir = normpath(joinpath(pwd(), dft.config["result_dir"]))
    if !isdir(result_dir)
        mkpath(result_dir)
    end

    df_dv = DataFrame(
                periods     = sort([k for k in dft.periods]),
                weights     = [dft.w[k] for k in sort([k for k in keys(dft.w)])],
                used_days   = [dft.u[k] for k in sort([k for k in keys(dft.u)])])
    writetable(joinpath(result_dir, "decision_variables.csv"), df_dv, separator=';')

    df_dv_s = deepcopy(df_dv[df_dv[:weights] .> 0, :])
    writetable(joinpath(result_dir, "decision_variables_short.csv"), df_dv_s, separator=';')

    ##################################################################################
    # Resulting time series
    ##################################################################################
    df_ts = DataFrame()
    period_idx = df_dv[:used_days] .== 1
    log("debug", "Period index of selected days: $period_idx")

    for ts in values(dft.time_series)
        df_ts[Symbol(ts.name)] =  ts.matrix_full[period_idx,:]'[:]
    end
    writetable(joinpath(result_dir, "resulting_profiles.csv"), df_ts, separator=';')

    ##################################################################################
    # Copy of config-file
    ##################################################################################
    cp(dft.config_file, joinpath(result_dir, "config_file.yaml"),remove_destination=true)
end

##################################################################################
# Create Plots of curves and reduced curves
##################################################################################
function createPlots(dft::DaysFinderTool)
    result_dir = normpath(joinpath(pwd(), dft.config["result_dir"], "plots"))
    if !isdir(result_dir)
        mkdir(result_dir)
    end

    w = [dft.w[k] for k in sort([k for k in keys(dft.w)])]
    u = [dft.u[k] for k in sort([k for k in keys(dft.w)])]

    periods_with_weight = [i for (i,k) in enumerate(sort([k for k in keys(dft.w)])) if dft.w[k]>0]
    weights = [ww for ww in w if ww > 0]

    for ts in values(dft.time_series)
        df = DataFrame(x=Float64[],y=Float64[], legend=String[])

        x = [x for x in range(1,length(ts.data))/length(ts.data)*100.]
        y = sort(ts.data, rev=true)
        append!(df, DataFrame(x=x, y=y, legend="original"))

        y = ts.matrix_full[periods_with_weight,:]'[:]
        x = (weights * ones(1, size(ts.matrix_full)[2]))'[:] / length(ts.data) * 100.
        df2 = sort(DataFrame(x=x, y=y, legend="reduced"), cols=:y, rev=true)
        df2[:x] = cumsum(df2[:x])
        append!(df, df2)


        p = plot(df, x=:x,y=:y,color=:legend, Geom.line,
                Guide.xlabel("Duration [-]"), Guide.ylabel("Curve"), Guide.title(ts.name),
                Coord.Cartesian(xmin=0,xmax=100))
        file_pdf = joinpath(result_dir, @sprintf "%s.pdf" ts.name)
        draw(PDF(file_pdf, 27cm, 21cm), p)
    end
end


##################################################################################
# Run Optimization
##################################################################################
function runDaysFinderToolDefault(dft::DaysFinderTool)
    @printf "%s\n" repeat("=", 100)
    @printf "%s\n" repeat("=", 100)

    dft.AREA_ERROR_TOLERANCE = Dict{String,Float64}()
    for c in dft.curves
        dft.AREA_ERROR_TOLERANCE[c] = dft.config["area_error_tolerance"]
    end

    u_val = Dict{String,Int}()
    w_val = Dict{String,Int}()
    for p in dft.periods
        u_val[p] = 0.
        w_val[p] = 0.
    end

    ##################################################################################
    # Model
    ##################################################################################
    # if dft.config["solver"]["mip"] == "CPLEX"
    #     m = Model(
    #         solver=CplexSolver(
    #                     CPX_PARAM_SCRIND=1,
    #                     CPX_PARAM_TILIM=dft.config["solver"]["ResLim"],
    #                     CPX_PARAM_THREADS=dft.config["solver"]["Threads"]))

    if dft.config["solver"]["mip"] == "GLPK"
        m = Model(solver=GLPKSolverMIP(presolve=true, msg_lev=GLPK.MSG_ON))
        # getSolverModel(m).params.tm_lim=15


    elseif dft.config["solver"]["mip"] == "CBC"
        m = Model(
                solver=CbcSolver(
                            seconds=dft.config["solver"]["ResLim"],
                            logLevel=1,
                            threads=32))

    elseif dft.config["solver"]["mip"] == "Gurobi"
        m = Model(
                solver=GurobiSolver(
                    Presolve=0,
                    TimeLimit=dft.config["solver"]["ResLim"],
                    Threads=dft.config["solver"]["Threads"]
                    )
                )

    else
        error()
        return
    end
    ##################################################################################
    # Variables
    ##################################################################################
    @variable(m, u[p in dft.periods], Bin, start=u_val[p])
    @variable(m, w[p in dft.periods] >= 0, start=w_val[p])

    @variable(m, area_error[c in dft.curves] >= 0)
    @variable(m, error[c in dft.curves, b in dft.bins] >= 0)


    # Mandatory Days
    m_d = String[]
    for ts in values(dft.time_series)
        mandatory_periods = get_mandatory_periods(ts, dft)
        append!(m_d, mandatory_periods)
        for p in mandatory_periods
            JuMP.fix(u[p],1)
            setlowerbound(w[p], 0.1)
        end
    end
    m_d = unique(m_d)

    # Starting Values
    p_start, w_start = define_random_start_point(dft, m_d)
    for p in p_start
        log("debug", "Set starting value for $p")
        log("debug", "Set weight to $(w_start[p])")
        if getcategory(u[p]) != :Fixed
            setvalue(u[p], 1)
        end
        setvalue(w[p], w_start[p])
    end

    ##################################################################################
    # Objective
    ##################################################################################
    @objective(m, Min, sum(dft.WEIGHT_DC[c] * sum(error[c,b] for b in dft.bins) for c in dft.curves))

    ##################################################################################
    # Constraints
    ##################################################################################
    # Defining error as absolute value
    @constraint(m, error_eq1[c in dft.curves, b in dft.bins],
        error[c,b] >= + dft.L[c,b] - sum( w[p] / dft.N_total_periods * dft.A[c,p,b] for p in dft.periods)
    )
    @constraint(m, error_eq2[c in dft.curves, b in dft.bins],
        error[c,b] >= - dft.L[c,b] + sum( w[p] / dft.N_total_periods * dft.A[c,p,b] for p in dft.periods)
    )

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
    log("debug","minimum weight set to 0.1")
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

    stat = solve(m)
    u_val = getvalue(u)
    w_val = getvalue(w)
    area_error_val = getvalue(area_error)


    if stat in [:Optimal, :UserLimit]
        u = getvalue(u)
        dft.u = Dict{String,Int}()
        for k in keys(u)
            dft.u[k[1]] = round(Int,abs(u[k[1]]))
        end

        w = getvalue(w)
        dft.w = Dict{String,Int}()
        for k in keys(w)
            dft.w[k[1]] = round(abs(w[k[1]]),4)
        end

    else
        @show stat
    end
    return stat
end