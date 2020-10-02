############################## Copyright (C) 2019  #############################
#       The content of this file is VITO (Vlaamse Instelling voor              #
#       Technologisch Onderzoek  N.V.) proprietary.                            #
################################################################################
module RepresentativeDaysFinders

    using Dates
    using TimeZones
    using DataFrames
    using CSV
    using YAML                              # -> read config-file
    using Combinatorics                     # -> functions for finding combinations
    using StatsBase                         # -> function to fit into Histogram bins
    using JuMP                              # -> optimization suite
    using Interpolations
    using Plots; gr()
    using JSON # output config file
    using SparseArrays # For the ordering variable
    using Cbc # For tests
    using GLPK # For tests
    using Gurobi # For tests
    using LinearAlgebra


    include("time_series/TimeSeries.jl")
    include("DaysFinderTools.jl")
    include("EmptyContainer.jl")
    include("time_series/TimeSeriesFunctions.jl")

    ###########################################################
    # functionality
    ###########################################################
    include("time_series/get_mandatory_periods.jl")
    include("output/create_plots.jl")
    include("output/writeOutResults.jl")
    include("util.jl")

    ###########################################################
    # data processing
    ###########################################################
    include("time_series/dataprocessing.jl")

    ############################################################
    # Default method to run tool
    ############################################################
    export findRepresentativeDays, ENTSOEcsv2dataframe, writeOutResults, DaysFinderTool, populateDaysFinderTool!, create_plots, makeDaysFinderToolModel, makeReOrderingDaysFinderTool, makeDCErrorOnlyDaysFinderToolModel, optimizeDaysFinderTool

    function findRepresentativeDays(
        dft::DaysFinderTool,
        optimizer_factory=GLPK.Optimizer
        )
        # Safety measure against having too many periods
        method = try_get_val(dft.config["solver"], "Method", "ordering")
        max_periods = 365
        if length(dft.periods) > max_periods && method == "ordering"
            @warn """
                Number of periods is greater than $max_periods and ordering solution method has been selected. This will likely crash Julia.
            """
        end

        # Make model
        if dft.config["solver"]["Method"] == "reorder"
            makeReOrderingDaysFinderTool(dft, optimizer_factory)
            stat = optimizeDaysFinderTool(dft)
        elseif dft.config["solver"]["Method"] == "squared reorder"
            makeReOrderingDaysFinderTool(dft, optimizer_factory)
            stat = optimizeDaysFinderTool(dft)
        elseif dft.config["solver"]["Method"] == "iterative_bins"
            binsVec = try_get_val(
                dft.config["solver"], "Bins", [10,20,40]
            )
            boundType = try_get_val(
                dft.config["solver"], "BoundType", "LowerBound"
            )
            for b in binsVec
                println("-"^80)
                println("-"^80)
                println("Number of bins used for DC error for this iteration: $b")
                println("-"^80)
                dft.config["number_bins"] = b
                if b != binsVec[1]
                    dft.config["hot_start_values"] = dft.v
                    if boundType == "ObjectiveLowerBound"
                        dft.config["objective_lower_bound"] =
                            objective_bound(dft.m)
                    elseif boundType == "ObjectiveValue"
                        dft.config["objective_lower_bound"] =
                            objective_value(dft.m)
                    end
                end
                populateDaysFinderTool!(dft)
                makeDaysFinderToolModel(dft, optimizer_factory)
                stat = optimizeDaysFinderTool(dft)
            end
        elseif dft.config["solver"]["Method"] == "DC_error_only"
            makeDCErrorOnlyDaysFinderToolModel(dft, optimizer_factory)
            fix_periods!(dft)
            stat = optimizeDaysFinderTool(dft)
        elseif dft.config["solver"]["Method"] == "ordering"
            makeDaysFinderToolModel(dft, optimizer_factory)
            fix_periods!(dft)
            stat = optimizeDaysFinderTool(dft)
        elseif dft.config["solver"]["Method"] == "squared error"
            makeSquaredErrorDaysFinderToolModel(dft, optimizer_factory)
            fix_periods!(dft)
            stat = optimizeDaysFinderTool(dft)
        elseif dft.config["solver"]["Method"] == "grow algorithm"
            @debug "Can't fix periods with this solution method"
            RP = solveGrowAlgorithmForPeriodSelection(dft, optimizer_factory)
            dft.rep_periods = RP
            makeReOrderingDaysFinderTool(dft, optimizer_factory)
            stat = optimizeDaysFinderTool(dft)
        elseif dft.config["solver"]["Method"] in ("chronological time period clustering", "hierarchical clustering")
            timePeriodClustering(dft)
            stat = MOI.OPTIMAL
        else
            error("Please specify solver method")
        end


        # write out results
        # TODO: The checks below could be a bit more robust...
        save_results = try_get_val(dft.config, "save_results", true)
        create_plots_bool = get(dft.config, "create_plots", 1)
        found_solution = (
            dft.config["solver"]["Method"] in (
                "chronological time period clustering",
                "hierarchical clustering"
                )
            ||
            has_values(dft.m)
        )
        if stat in [MOI.OPTIMAL, MOI.TIME_LIMIT] && found_solution && save_results
            writeOutResults(dft)
            if create_plots_bool in (1,true)
                create_plots(dft)
            end
        end
        return dft
    end

    function findRepresentativeDays(config_file::String, optimizer_factory)
        @info("Start RepresentativeDaysFinder")
        if isfile(config_file)
            dft = DaysFinderTool(config_file, populate_entries=true)
            findRepresentativeDays(dft, optimizer_factory)
        else
            @error("Config-file not found $config_file")
        end
    end

end
