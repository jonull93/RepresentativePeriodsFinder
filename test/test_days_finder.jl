using RepresentativePeriodsFinder
RPF = RepresentativePeriodsFinder
using JuMP
using Cbc

# Create PeriodsFinder
config_file = normpath(joinpath(@__DIR__, "input_data", "default.yaml"))
pf = PeriodsFinder(config_file, populate_entries=true)

# Delete entries for time series and duration curve error
delete!(pf.config["method"]["optimization"], "time_series_error")
delete!(pf.config["method"]["options"], "ordering_error")
delete!(pf.config["method"], "clustering")

# Make life a bit easier for the optimiser, reduce the number of days
pf.config["method"]["options"]["total_periods"] = 20

# Setup optimizer
opt = optimizer_with_attributes(Cbc.Optimizer, "seconds" => 60)
# using Gurobi
# opt = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => 60)

# Create and solve optimization problem
# find_representative_periods(pf, optimizer=opt, reset=true)

# Try absolute error (this is the only possibility with Cbc)
pf.config["method"]["optimization"]["duration_curve_error"]["type"] = "absolute"
find_representative_periods(pf, optimizer=opt, reset=true)