using RepresentativeDaysFinders
using JuMP
using Gurobi
timeLimitVal = 60*1 # in seconds
config_file = normpath(joinpath(@__DIR__, "Elia_2017.yaml"))

dft = DaysFinderTool(config_file)
dft.config["number_days"] = 100 # Number of periods to choose
dft.config["number_days_total"] = 8760 # Number of periods in year
dft.config["timesteps_per_period"] = 1 # Number of timesteps per period
dft.config["result_dir"] = joinpath(@__DIR__, "results_rep_hours_finder") # Where to save results - this directory is relative to the .yaml file
# Specifying an absolute directory (e.g C://Users/me/Desktop/results) should also work, but won't work on Mac or Linux (if you're using Windows).
dft.config["solver"]["Method"] = "DC_error_only" # This is the original rep days finder methodology
dft.config["integral_weights"]=false

# Populate the days finder tool. This adds the time series to it.
populateDaysFinderTool!(dft)

# Howto fix #num_max_periods in your data:
num_max_periods = 30
max_ix = sortperm(reshape(dft.time_series["Total_load"].matrix_full, size(dft.time_series["Total_load"].matrix_full)[1]), rev=true)
dft.config["fixed_periods"] = max_ix[1:num_max_periods]


# Make the model first
optimizer_factory = RepresentativeDaysFinders.optimizer_with_attributes(
    Gurobi.Optimizer, "TimeLimit" => timeLimitVal
)
RepresentativeDaysFinders.makeDCErrorOnlyDaysFinderToolModel(
    dft, optimizer_factory
)

# alternative you can also fix weights as
max_ix = sortperm(reshape(dft.time_series["Total_load"].matrix_full, size(dft.time_series["Total_load"].matrix_full)[1]), rev=true)[1:num_max_periods]
u = dft.misc[:u]
w = dft.misc[:w]
for ix in max_ix
    fix(u[ix], 1, force=true)
    fix(w[ix], 1.0, force=true)
end

# Optimise here
stat = optimizeDaysFinderTool(dft)

# Write out the results
writeOutResults(dft)

# Create plots (if not done already)
create_plots(dft)
