# # Selecting representative hours using optimisation
# In this example we'll see how to select representative hours from a year using the optimisation based approach.

# First we create an instance of a PeriodsFinder. We'll use the default `.yaml` configuration file and modify it, but you could just edit it directly.
using RepresentativePeriodsFinder, PrettyPrint, JuMP;
RPF = RepresentativePeriodsFinder;
config_file = normpath(joinpath(@__DIR__, "input_data", "default.yaml"));
pf = PeriodsFinder(config_file, populate_entries=false);
pprintln(pf.config)

# Delete entries for clustering, time series and ordering error. This way we'll select hours.
delete!(pf.config["method"]["optimization"], "time_series_error");
delete!(pf.config["method"]["options"], "ordering_error");
# delete!(pf.config["method"]["optimization"], "duration_curve_error");
delete!(pf.config["method"], "clustering");

# temporary deletion of all timeseries except for solar
for ts_name in ["LFW", "Load", "Residual Load"] 
    delete!(pf.config["time_series"], ts_name)
end

# Use absolute error (this is the only possibility with MILP solvers)
pf.config["method"]["optimization"]["duration_curve_error"]["type"] = "squared";
pf.config["method"]["optimization"]["duration_curve_error"]["number_bins"] = 1000;
# pf.config["method"]["optimization"]["time_series_error"]["type"] = "absolute";

# Choose the number of representative hours to select
n_rep = 23;
pf.config["method"]["options"]["representative_periods"] = n_rep;

# Make life a bit easier for the optimiser, let's reduce the number of hours from which we're selecting, i.e. we're shortening the year to the first couple of days
pf.config["method"]["options"]["total_periods"] = 24;

# Change period length to 1 (24 = day selection)
pf.config["method"]["options"]["time_steps_per_period"] = 1;

# Turn off automatic plot creation
pf.config["results"]["create_plots"] = false;

# Give a very high weight to the solar time series
pf.config["time_series"]["LFS"]["weight"] = 100

# Change the result directory name.
script_name = splitext(splitdir(@__FILE__)[2])[1];
result_dir = joinpath(@__DIR__, "results", script_name);
pf.config["results"]["result_dir"] = result_dir;

# This is what the resulting configuration file looks like now:
pprintln(pf.config);

# Populate the entries now (could have also done this initially)
populate_entries!(pf)

# display(Plots.plot(RPF.get_discretised_duration_curve(pf, "LFS")))
# y = RPF.get_bin_interval_values(pf, "LFS")
# y = [mean([y[i+1], y[i]]) for i in eachindex(y)[1:end-1]]
# y = RPF.get_set_of_bins(pf)
# L = RPF.get_discretised_duration_curve(pf, "LFS")
# display(Plots.plot(L))
# nts = RPF.get_normalised_time_series_values(pf, "LFS")

# Setup optimizer. Let's use HiGHS for a change.
# opt = optimizer_with_attributes(HiGHS.Optimizer, "time_limit" => 60.0, "output_flag" => true);
opt = optimizer_with_attributes(Gurobi.Optimizer, "time_limit" => 60.0, "output_flag" => true);

# Find the representative hours
find_representative_periods(pf, optimizer=opt, reset=true);

# Plot the duration curves
for (ts_name, ts) in RPF.get_set_of_time_series(pf)
    display(create_duration_curve(pf, ts_name))
end

# create_synthetic_time_series_plot(pf, "LFS")