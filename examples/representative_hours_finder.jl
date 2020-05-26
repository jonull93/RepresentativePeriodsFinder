using RepresentativeDaysFinders
using JuMP
using Gurobi
timeLimitVal = 60 # in seconds
config_file = normpath(joinpath(@__DIR__, "input_data", "Elia_2017.yaml"))

dft = DaysFinderTool(config_file)
dft.config["number_days"] = 100 # Number of periods to choose
dft.config["number_days_total"] = 8760 # Number of periods in year
dft.config["timesteps_per_period"] = 1 # Number of timesteps per period
dft.config["result_dir"] = joinpath("..", "results", "rep_hours_finder") # Where to save results - this directory is relative to the .yaml file
# Specifying an absolute directory (e.g C://Users/me/Desktop/results) should also work, but won't work on Mac or Linux (if you're using Windows).
dft.config["solver"]["Method"] = "DC_error_only" # This is the original rep days finder methodology

# Populate the days finder tool. This adds the time series to it.
populateDaysFinderTool!(dft)

# Make the model first
optimizer_factory = optimizer_with_attributes(
    Gurobi.Optimizer, "TimeLimit" => timeLimitVal
)
RepresentativeDaysFinders.makeDCErrorOnlyDaysFinderToolModel(
    dft, optimizer_factory
)

# You can fix any hours / days you want here
# As an example, here I fix hour 1
u = dft.misc[:u]
fix(u[1], 1, force=true)

# Optimise here
stat = optimizeDaysFinderTool(dft)

# Write out the results
writeOutResults(dft)

# Create plots (if not done already)
create_plots(dft)
