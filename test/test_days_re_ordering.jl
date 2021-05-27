# # Re-ordering days in a year
# This script illustrates how to re-order days in a year, for example in order to go from 10 representative days to a continuous year of 365 (representative) days.

# First we create an instance of a PeriodsFinder. Here we are using the default `.yaml` configuration file, but keep in mind that a lot of the following lines of code could be avoided by creating your own custom file.
using RepresentativePeriodsFinder
RPF = RepresentativePeriodsFinder
config_file = normpath(joinpath(@__DIR__, "input_data", "default.yaml"))
pf = PeriodsFinder(config_file, populate_entries=true)

# Change the result directory name
script_name = splitext(splitdir(@__FILE__)[2])[1]
result_dir = joinpath(@__DIR__, "results", script_name)
pf.config["results"]["result_dir"] = result_dir

# We will only concern ourselves with mapping load, just to make things simple.
delete!(pf.config["time_series"], "LFW")
delete!(pf.config["time_series"], "LFS",)
delete!(pf.config["time_series"], "Residual Load")

# Next we define the total number of periods and the number of representative periods. This script was run without access to a professional solver such as CPLEX or Gurobi, so we will limit ourselves to creating a "year" of 20 days. Further on we will use an alternative formulation which will be able to create a full year of 365 days.
N_total = 20
N_rep = 10
pf.config["method"]["options"]["total_periods"] = 20
pf.config["method"]["options"]["representative_periods"] = 10

# We define our representative periods, where each element in the vector is a day in the year.
rep_periods = [1, 3, 4, 5, 7, 12, 15, 16, 18, 19]

# We populate the results of the selection variable, `pf.u`, with our representative days.
pf.u = zeros(Bool, N_total)
pf.u[rep_periods] .= 1

# The following lines define the selection and ordering method, specifically we will use an optimisation method which minimises only the squared (as opposed to absolute) time series error. 
delete!(pf.config["method"], "clustering")
delete!(pf.config["method"]["optimization"], "duration_curve_error")
pf.config["method"]["optimization"]["time_series_error"]["type"] = "absolute"

# We set binary ordering to false, which means that a non-representative day can be represented by a linear combination of representative days.
pf.config["method"]["optimization"]["binary_ordering"] = false

# We set the mandatory periods to be selected to our representative periods
pf.config["method"]["options"]["mandatory_periods"] = rep_periods

# Setup the optimizer and solve to re-order the days.
using Ipopt
opt = optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 100)
find_representative_periods(pf, optimizer=opt, reset=false)

# We can plot the resulting synthetic time series. Note how days 1 and 3, which are representative, match up exactly while other days must be approximated. 
using Dates
timestamps = Dict("Load" => DateTime(1970,1,1):Hour(1):DateTime(1970,1,7))
create_synthetic_time_series_plots(pf; timestamps=timestamps)
# ![synth_load_time_series](joinpath(result_dir, "Load_synthetic_time_series.svg))

# We can also write out the this synthetic time series to `pf.config["results"]["result_dir"].
write_out_synthetic_timeseries(pf, timestamps=timestamps)