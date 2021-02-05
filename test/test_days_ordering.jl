using RepresentativePeriodsFinder
RPF = RepresentativePeriodsFinder
using JuMP
using Cbc

# Create PeriodsFinder
config_file = normpath(joinpath(@__DIR__, "input_data", "default.yaml"))
pf = PeriodsFinder(config_file, populate_entries=true)

# Delete entries for time series error and ordering error
delete!(pf.config["method"]["optimization"], "time_series_error")
delete!(pf.config["method"]["optimization"], "duration_curve_error")

# Can edit the number of representative days to be chosen here
pf.config["method"]["options"]["representative_periods"] = 64

# Make optimisation model
m = RPF.make_periods_finder_model!(pf)

# Optimise it
RPF.optimize_periods_finder_model!(pf, m)

# Create the plots
RPF.create_plots(pf)

# Delete results
rm(RPF.get_abspath_to_result_dir(pf), recursive=true)