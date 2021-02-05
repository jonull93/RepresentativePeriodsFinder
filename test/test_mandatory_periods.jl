using RepresentativePeriodsFinder
RPF = RepresentativePeriodsFinder
using JuMP
using Cbc
using StatsBase
# using Gurobi

# Create PeriodsFinder
config_file = normpath(joinpath(@__DIR__, "input_data", "default.yaml"))
pf = PeriodsFinder(config_file, populate_entries=true)

# Delete clustering entries and ensure problem is linear
delete!(pf.config["method"], "clustering")
delete!(pf.config["method"]["optimization"], "time_series_error")
delete!(pf.config["method"]["optimization"], "duration_curve_error")
pf.config["method"]["optimization"]["binary_ordering"] = true

# Setup optimizer
opt = optimizer_with_attributes(Cbc.Optimizer)

# Special case of all periods set to mandatory
n_rep = RPF.get_number_of_representative_periods(pf)
pf.config["method"]["options"]["mandatory_periods"] = [i for i in 1:n_rep]
find_representative_periods(pf, optimizer=opt)
@assert all(pf.u[1:n_rep] .== 1) 

# Only on the first day of the year
pf.config["method"]["options"]["mandatory_periods"] = [1]
find_representative_periods(pf, optimizer=opt)
@assert pf.u[1] == 1

# Same for clustering
delete!(pf.config["method"], "optimization")
pf.config["method"]["clustering"] = Dict(
    "type" => "hierarchical"
)
find_representative_periods(pf)
@assert pf.u[1] == 1

# Clean up
rm(RPF.get_abspath_to_result_dir(pf), recursive=true)