# using RepresentativePeriodsFinder
# RPF = RepresentativePeriodsFinder

# # Create PeriodsFinder
# config_file = normpath(joinpath(@__DIR__, "input_data", "default.yaml"))
# pf = PeriodsFinder(config_file, populate_entries=true)

# # Delete entries for time series and duration curve error
# delete!(pf.config["method"], "optimization")

# # Create and solve optimization problem
# find_representative_periods(pf, reset=true)

# # Clean up
# rm(RPF.get_abspath_to_result_dir(pf), recursive=true)