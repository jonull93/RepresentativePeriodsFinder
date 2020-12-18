using RepresentativePeriodsFinder
using JuMP
using Cbc

timeLimitVal = 60 # in seconds
config_file = normpath(joinpath(@__DIR__, "input_data", "Elia_2017.yaml"))

# Change the value of mult below to go from days to half days, ... down to hours. Note that it's REALLY slow for hours... and that's if nothing crashes.
mult = 8
timeLimitVal = timeLimitVal*mult^2 # Just to be on the safe side

dft = PeriodsFinder(config_file)
dft.config["number_days"] = 8*mult # Number of periods to choose
dft.config["number_days_total"] = 365*mult # Number of periods in year
dft.config["timesteps_per_period"] = Int(24//mult) # Number of timesteps per period
dft.config["duration_curve_error_weight"] = 0.5 # weight to your choice
dft.config["time_series_error_weight"] = 0.5 # same here
dft.config["time_series_error_matrix_type"] = "absolute" # Not too important
dft.config["result_dir"] = joinpath("..", "results", "rep_days_order") # Where to save results - this directory is relative to the .yaml file
# Specifying an absolute directory (e.g C://Users/me/Desktop/results) should also work, but won't work on Mac or Linux (if you're using Windows).

# Setting below tells dft to find and order representative periods using an integer problem, not the original representative days finder.
dft.config["solver"]["Method"] = "ordering"

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
