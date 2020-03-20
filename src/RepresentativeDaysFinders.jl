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
    export findRepresentativeDays, ENTSOEcsv2dataframe, writeOutResults, DaysFinderTool, populateDaysFinderTool!, create_plots

    function findRepresentativeDays(dft::DaysFinderTool, optimizer_factory)
        # Optimize
        if dft.config["solver"]["Method"] == "BruteForce"
            stat = runDaysFinderToolBruteForce(dft, optimizer_factory)
        elseif dft.config["solver"]["Method"] == "Iterative"
            stat = runDaysFinderToolIterativeBounderaries(dft, optimizer_factory)
        else
            stat = runDaysFinderToolDefault(dft, optimizer_factory)
        end

        # write out results
        if stat in [MOI.OPTIMAL, MOI.TIME_LIMIT]
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
