using RepresentativePeriodsFinder
RPF = RepresentativePeriodsFinder
using JuMP
using Ipopt
using StatsBase
# using Gurobi

# Create PeriodsFinder
config_file = normpath(joinpath(@__DIR__, "input_data", "default.yaml"))
pf = PeriodsFinder(config_file, populate_entries=true)

# Select random days as representative
pf.config["method"]["options"]["representative_periods"] = 10
N_total = RPF.get_number_of_periods(pf)
N_rep = RPF.get_number_of_representative_periods(pf)
pf.u = zeros(Bool, N_total)
rep_periods = sample(1:365, N_rep, replace=false)
pf.u[rep_periods] .= 1

# Delete clustering entries and ensure problem is linear
delete!(pf.config["method"], "clustering")
delete!(pf.config["method"]["optimization"], "duration_curve_error")
pf.config["method"]["optimization"]["binary_ordering"] = false
pf.config["method"]["options"]["mandatory_periods"] = rep_periods

# Setup optimizer
opt = optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 100)
# opt = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => 60)

# Create and solve optimization problem - with reset=false for re-ordering!
find_representative_periods(pf, optimizer=opt, reset=false)

# Try absolute error
pf.config["method"]["optimization"]["time_series_error"]["type"] = "absolute"
find_representative_periods(pf, optimizer=opt, reset=false)

# Try writing these out to the results
timestamps = Dict("Load" => DateTime(1970,1,1):Hour(1):DateTime(1970,1,2))
create_synthetic_time_series_plots(pf, timestamps=timestamps)
write_out_synthetic_timeseries(pf, timestamps=timestamps)

# Clean up
rm(RPF.get_abspath_to_result_dir(pf), recursive=true)