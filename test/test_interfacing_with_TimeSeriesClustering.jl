# # Selecting representative hours using optimisation

# In this example we'll see convert a `PeriodsFinder` to a `TimeSerisClustering.ClustData` type, so that we can leverage the clustering algorithms from the [TimeSeriesClustering.jl package](https://github.com/holgerteichgraeber/TimeSeriesClustering.jl).

# Keep in mind that you will probably have to use a fork of that package for now due to package compatibility issues:

using Pkg
pkg"add https://github.com/junglegobs/TimeSeriesClustering.jl.git"

# Make the periods finder
using RepresentativePeriodsFinder, TimeSeriesClustering, TimeSeries
RPF = RepresentativePeriodsFinder;
config_file = RPF.datadir("default.yaml");
pf = PeriodsFinder(config_file; populate_entries=true);

# Let's define the following functions to convert between these two packages:

function Base.convert(::Type{TimeSeriesClustering.ClustData}, pf::PeriodsFinder)
    region = "none"
    years = [2016]
    K = RPF.get_number_of_periods(pf)
    T = RPF.get_number_of_time_steps_per_period(pf)
    data = Dict(
        ts_name => reshape(values(ts)[1:(T * K)], T, K) for
        (ts_name, ts) in RPF.get_set_of_time_series(pf)
    )
    weights = [1.0 for i in 1:K]
    mean = Dict(
        ts_name => [0.0 for t in 1:T] for
        ts_name in RPF.get_set_of_time_series_names(pf)
    )
    sdv = Dict(
        ts_name => [1.0 for t in 1:T] for
        ts_name in RPF.get_set_of_time_series_names(pf)
    )
    delta_t = [1.0 for t in 1:T, i in 1:K]
    k_ids = [i for i in 1:K]
    return ClustData(
        region, years, K, T, data, weights, mean, sdv, delta_t, k_ids
    )
end

# Let's convert this to a `ClustData` type:
cd = convert(TimeSeriesClustering.ClustData, pf);

# Let's run the hierarchical clustering (and time it as well, to see how it compares with the inbuilt clustering algorithm):
t_1 = @elapsed cr = run_clust(
    cd; method="hierarchical", representation="centroid", n_clust=8, n_init=1
)

# Converting back to a `PeriodsFinder` type is less straightforward, as the above clustering does not save which days where selected i.e. there is no `u` variable. We would compare the values of a cluster to those in each period in the original time series to find out which is which.

# In any case, it may be interesting to see if the above clustering was much faster than the inbuilt `RepresentativePeriodsFinder` clustering algorithm:
delete!(pf.config["method"], "optimization")
delete!(pf.config["method"]["options"]["ordering_error"], "ord_err_2")
pf.inputs[:ordering_error_functions]["ord_err_1"] = (
   (x,y) -> (sum((x[i] - y[i])^2 for i in eachindex(x)))
)
t_2 = @elapsed find_representative_periods(pf)

# The answer is manifestly yes, even accounting for Julia's compilation times, which is good to know.

#jl Pkg.rm("TimeSeriesClustering")