############################## Copyright (C) 2019  #############################
#       The content of this file is VITO (Vlaamse Instelling voor              #
#       Technologisch Onderzoek  N.V.) proprietary.                            #
################################################################################
module RepresentativePeriodsFinder

    using CSV                               # -> read csv files
    using DataFrames                        # -> csv read into DataFrames
    using TimeZones                         # -> time zone conversion
    using Dates                             # -> date time parsing
    using TimeSeries                        # -> time series type
    using YAML                              # -> read config-file
    using Combinatorics                     # -> functions for finding combinations
    using StatsBase                         # -> function to fit into Histogram bins
    using JuMP                              # -> optimization suite
    using Interpolations                    # -> interpolating time series
    using Plots; gr()                       # -> output plots
    using JSON                              # -> output config file
    using SparseArrays                      # -> for the (sparse) ordering variable
    using Cbc                               # -> free solver for tests
    using LinearAlgebra                     # -> for e.g. transpose()
    using NamedTupleTools                   # -> for easy creation of NamedTuples

    # # PeriodsFinder type
    include("PeriodsFinder.jl")

    # # Input and data processing
    include("input/read_time_series.jl")
    include("input/process_time_series.jl")
    include("input/load_periods_finder.jl")

    # # Period finding methods
    # include("methods/optimisation")
    # include("methods/clustering")

    # # Outputs and results
    # include("output/create_plots.jl")
    # include("output/write_out_results.jl")

    # # Utility
    include("util/EmptyContainer.jl")
    include("util/util.jl")

    # Exported methods
    export find_representative_periods, ENTSOEcsv2dataframe, writeOutResults, PeriodsFinder, populate_days_finder!, create_plots, makePeriodsFinderModel, makeReOrderingPeriodsFinder, makeDCErrorOnlyPeriodsFinderModel, optimizePeriodsFinder

    # function find_representative_periods(
    #     dft::PeriodsFinder,
    #     optimizer_factory=Cbc.Optimizer
    #     )
    #     # Safety measure against having too many periods
    #     method = try_get_val(dft.config["solver"], "Method", "ordering")
    #     max_periods = 365
    #     if length(dft.periods) > max_periods && method == "ordering"
    #         @warn """
    #             Number of periods is greater than $max_periods and ordering solution method has been selected. This will likely crash Julia.
    #         """
    #     end

    #     # Make model
    #     if dft.config["solver"]["Method"] == "reorder"
    #         makeReOrderingPeriodsFinder(dft, optimizer_factory)
    #         stat = optimizePeriodsFinder(dft)
    #     elseif dft.config["solver"]["Method"] == "squared reorder"
    #         makeReOrderingPeriodsFinder(dft, optimizer_factory)
    #         stat = optimizePeriodsFinder(dft)
    #     elseif dft.config["solver"]["Method"] == "iterative_bins"
    #         binsVec = try_get_val(
    #             dft.config["solver"], "Bins", [10,20,40]
    #         )
    #         boundType = try_get_val(
    #             dft.config["solver"], "BoundType", "LowerBound"
    #         )
    #         for b in binsVec
    #             println("-"^80)
    #             println("-"^80)
    #             println("Number of bins used for DC error for this iteration: $b")
    #             println("-"^80)
    #             dft.config["number_bins"] = b
    #             if b != binsVec[1]
    #                 dft.config["hot_start_values"] = dft.v
    #                 if boundType == "ObjectiveLowerBound"
    #                     dft.config["objective_lower_bound"] =
    #                         objective_bound(dft.m)
    #                 elseif boundType == "ObjectiveValue"
    #                     dft.config["objective_lower_bound"] =
    #                         objective_value(dft.m)
    #                 end
    #             end
    #             populate_days_finder!(dft)
    #             makePeriodsFinderModel(dft, optimizer_factory)
    #             stat = optimizePeriodsFinder(dft)
    #         end
    #     elseif dft.config["solver"]["Method"] == "DC_error_only"
    #         makeDCErrorOnlyPeriodsFinderModel(dft, optimizer_factory)
    #         fix_periods!(dft)
    #         stat = optimizePeriodsFinder(dft)
    #     elseif dft.config["solver"]["Method"] == "ordering"
    #         makePeriodsFinderModel(dft, optimizer_factory)
    #         fix_periods!(dft)
    #         stat = optimizePeriodsFinder(dft)
    #     elseif dft.config["solver"]["Method"] == "squared error"
    #         makeSquaredErrorPeriodsFinderModel(dft, optimizer_factory)
    #         fix_periods!(dft)
    #         stat = optimizePeriodsFinder(dft)
    #     elseif dft.config["solver"]["Method"] == "grow algorithm"
    #         @debug "Can't fix periods with this solution method"
    #         RP = solveGrowAlgorithmForPeriodSelection(dft, optimizer_factory)
    #         dft.rep_periods = RP
    #         makeReOrderingPeriodsFinder(dft, optimizer_factory)
    #         stat = optimizePeriodsFinder(dft)
    #     elseif dft.config["solver"]["Method"] in ("chronological time period clustering", "hierarchical clustering")
    #         timePeriodClustering(dft)
    #         stat = MOI.OPTIMAL
    #     else
    #         error("Please specify solver method")
    #     end


    #     # write out results
    #     # TODO: The checks below could be a bit more robust...
    #     save_results = try_get_val(dft.config, "save_results", true)
    #     create_plots_bool = get(dft.config, "create_plots", 1)
    #     found_solution = (
    #         dft.config["solver"]["Method"] in (
    #             "chronological time period clustering",
    #             "hierarchical clustering"
    #             )
    #         ||
    #         has_values(dft.m)
    #     )
    #     if stat in [MOI.OPTIMAL, MOI.TIME_LIMIT] && found_solution && save_results
    #         writeOutResults(dft)
    #         if create_plots_bool in (1,true)
    #             create_plots(dft)
    #         end
    #     end
    #     return dft
    # end

    # function find_representative_periods(config_file::String, optimizer_factory)
    #     @info("Start RepresentativeDaysFinder")
    #     if isfile(config_file)
    #         dft = PeriodsFinder(config_file, populate_entries=true)
    #         find_representative_periods(dft, optimizer_factory)
    #     else
    #         @error("Config-file not found $config_file")
    #     end
    # end

end


