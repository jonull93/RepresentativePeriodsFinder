using RepresentativePeriodsFinder
using JuMP
using Cbc

# Create PeriodsFinder
config_file = normpath(joinpath(@__DIR__, "input_data", "default.yaml"))
pf = PeriodsFinder(config_file, populate_entries=true)

# Delete entries for time series error
delete!(pf.config["method"]["optimization"], "time_series_error")
delete!(pf.config["method"]["optimization"], "duration_curve_error")

# Make optimisation model
m = RepresentativePeriodsFinder.make_periods_finder_model!(pf)
