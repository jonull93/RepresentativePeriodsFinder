##################################################################################
# Author:   Hanspeter Höschle
# Date:     15/06/2017
##################################################################################

module RepresentativeDaysFinder

    # Method for logging tool output -> log("info", "MESSAGE")
    using Lumberjack.log

    ##################################################################################
    # type declaration for TimeSeries
    ##################################################################################
    using DataFrames
    include("TimeSeries.jl")

    ##################################################################################
    # type declaration and functions for DaysFinderTool
    ##################################################################################
    using YAML                              # -> read config-file

    using Combinatorics                     # -> functions for finding combinations

    using JuMP                              # -> optimization suite
    # using CPLEX                             # -> Solver interface: CPLEX
    using GLPK,GLPKMathProgInterface        # -> Solver interface: GLPK
    using Cbc                               # -> Solver interface: Cbc
    using Gurobi                            # -> Solver interface: Gurobi


    using Gadfly
    include("DaysFinderTools.jl")

    ##################################################################################
    # Functions for TimeSeries depending on type DaysFinderTool
    ##################################################################################
    include("TimeSeriesFunctions.jl")

    ##################################################################################
    # Default method to run tool
    ##################################################################################
    export findRepresentativeDays

    function findRepresentativeDays(config_file::String)
        log("info", "Start RepresentativeDaysFinder")
        if isfile(config_file)
            dft = DaysFinderTool(config_file)
            if dft.config["solver"]["Method"] == "BruteForce"
                stat = runDaysFinderToolBruteForce(dft)
            elseif dft.config["solver"]["Method"] == "Iterative"
                stat = runDaysFinderToolIterativeBounderaries(dft)
            else
                stat = runDaysFinderToolDefault(dft)
            end

            if stat in [:Optimal, :UserLimit]
                writeOutResults(dft)
                if dft.config["create_plots"] == 1
                    createPlots(dft)
                end
            end
        else
            log("error", "Config-file not found $config_file")
        end
    end

end
