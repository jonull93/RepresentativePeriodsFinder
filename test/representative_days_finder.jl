using RepresentativePeriodsFinder
using JuMP
using Cbc
timeLimitVal = 60 # in seconds
config_file = normpath(joinpath(@__DIR__, "input_data", "Elia_2017.yaml"))

dft = PeriodsFinder(config_file)
dft.config["number_days"] = 8 # Number of periods to choose
dft.config["number_days_total"] = 365 # Number of periods in year
dft.config["timesteps_per_period"] = 24 # Number of timesteps per period
dft.config["result_dir"] = joinpath("..", "results", "rep_hours_finder") # Where to save results - this directory is relative to the .yaml file
# Specifying an absolute directory (e.g C://Users/me/Desktop/results) should also work, but won't work on Mac or Linux (if you're using Windows).
dft.config["solver"]["Method"] = "DC_error_only" # This is the original rep days finder methodology

# Populate the days finder tool. This adds the time series to it.
populate_entries!(dft)

# Now we can find the representative days
find_representative_periods(dft,
    RepresentativePeriodsFinder.optimizer_with_attributes(
        Cbc.Optimizer, "TimeLimit" => timeLimitVal
    )
)

# Write out the results
writeOutResults(dft)

# Create plots (if not done already)
create_plots(dft)
