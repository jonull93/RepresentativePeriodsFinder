# # Ramping and correlation time series

# In the original paper by [Poncelet et al.](https://www.mech.kuleuven.be/en/tme/research/energy_environment/Pdf/wp-2015-10b.pdf) there was the possibility of using a correlation time series, e.g. for electrical load and solar power (which are highly correlated over the course of a day). Equally, it may be interesting to select representative periods based on a ramping or mileage time series (e.g. to select days with a high rate of change in the load for example).

# Eventually, this could be added as a core feature within `RepresentativePeriodsFinder.jl`, and indeed there are still two open issues regarding this (see [here](https://gitlab.kuleuven.be/UCM/representativeperiodsfinder.jl/-/issues/16) and [here](https://gitlab.kuleuven.be/UCM/representativeperiodsfinder.jl/-/issues/17)). This is not a priority however and will not be done unless someone (that could be you!) opens a pull request. Anyway, enough faffing about - let's get cracking!

# ## Ramping time series

# Let's create a `PeriodsFinder` instance:
using RepresentativePeriodsFinder, Plots;
RPF = RepresentativePeriodsFinder;
config_file = RPF.datadir("default.yaml");
pf = PeriodsFinder(config_file; populate_entries=false);

# Delete some time series to speed up the process
ts_orig = copy(pf.config["time_series"]);
for ts_name in ["LFW", "Load", "LFS"]
    delete!(pf.config["time_series"], ts_name)
end
populate_entries!(pf);

# Let's now create a ramping time series based on the residual load time series. Using `deepcopy` instead of `copy` is extremely important if you want to avoid weird bugs!
pf.time_series["RL Ramp"] = diff(deepcopy(pf.time_series["Residual Load"]));

# The above time series is 8783 entries long, but that's not as an issue - the configuration file will only select the first 8,760 entries. However, it starts an hour earlier than "Residual Load" - this is an issue! Let's fix it.
timestamp(pf.time_series["RL Ramp"]) .-= Hour(1)

# We also need to change some settings in the "RL Ramp" meta dictionary:
meta(pf.time_series["RL Ramp"])["name"] = "RL Ramp";
meta(pf.time_series["RL Ramp"])["weight"] = 0.00001;
meta(pf.time_series["Residual Load"])

# TODO: allow for weights in clustering algorithm! Also the weights aren't working for some reason?

# And also add an entry for it in the config dictionary. This avoids some unfortunate bugs.
pf.config["time_series"]["RL Ramp"] = Dict()

delete!(pf.config["method"]["optimization"], "time_series_error")
delete!(pf.config["method"]["options"], "ordering_error")
delete!(pf.config["method"], "clustering")
pf.config["method"]["optimization"]["duration_curve_error"]["type"] = "absolute"

# # Make this a clustering based selection:
# delete!(pf.config["method"], "optimization")
# delete!(pf.config["method"]["options"]["ordering_error"], "ord_err_2"); # remove from set
# pf.inputs[:ordering_error_functions]["ord_err_1"] = (
#     (x, y) -> (sum((x[i] - y[i])^2 for i in eachindex(x)))
# );

# And now let's find us some representative days:
find_representative_periods(
    pf; optimizer=optimizer_with_attributes(HiGHS.Optimizer), reset=false
)
# find_representative_periods(pf; reset=false)

# And what do these duration curves look like then?
plt_vec = [
    create_duration_curve(pf, ts_name; line=:solid, marker=:none) for
    (ts_name, ts) in RPF.get_set_of_time_series(pf)
]
display(Plots.plot(plt_vec...; layout=(1, 2)))