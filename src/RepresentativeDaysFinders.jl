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
    export findRepresentativeDays, ENTSOEcsv2dataframe, writeOutResults, DaysFinderTool, populateDaysFinderTool!, create_plots, makeDaysFinderToolModel, postOptimizationOfOrdering

    function findRepresentativeDays(dft::DaysFinderTool, optimizer_factory)
        # Safety measure against having too many periods
        method = try_get_val(dft.config["solver"], "Method", "ordering")
        max_periods = 500
        if length(dft.periods) > max_periods && method == "ordering"
            error("Number of periods is greater than $max_periods - this will probably crash Julia.")
        end

        # Make model
        if dft.config["solver"]["Method"] == "reorder"
            try
                makeReOrderingDaysFinderTool(dft, optimizer_factory)
                stat = optimizeDaysFinderTool(dft)
            catch
                error("Could not reorder periods. Did you forget to solve to find periods?")
            end
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
            stat = optimizeDaysFinderTool(dft)
        elseif dft.config["solver"]["Method"] == "ordering"
            makeDaysFinderToolModel(dft, optimizer_factory)
            stat = optimizeDaysFinderTool(dft)
        else
            error("Please specify solver method")
        end


        # write out results
        save_results = try_get_val(dft.config, "save_results", true)
        if stat in [MOI.OPTIMAL, MOI.TIME_LIMIT] && save_results
            writeOutResults(dft)
            if dft.config["create_plots"] == 1
                create_plots(dft)
            end
        end
        return dft
    end

    function findRepresentativeDays(config_file::String, optimizer_factory)
        @info("Start RepresentativeDaysFinder")
        if isfile(config_file)
            dft = DaysFinderTool(config_file)
            findRepresentativeDays(dft, optimizer_factory)
        else
            @error("Config-file not found $config_file")
        end
    end

end
