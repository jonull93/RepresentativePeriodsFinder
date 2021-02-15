using RepresentativePeriodsFinder, Dates, Cbc
RPF = RepresentativePeriodsFinder

# Create PeriodsFinder
config_file = normpath(joinpath(@__DIR__, "input_data", "default.yaml"))
pf = PeriodsFinder(config_file, populate_entries=true)

# Delete entries for optimization, add intermediate periods
# and use only 1 ordering error
delete!(pf.config["method"], "optimization")
pf.config["method"]["clustering"]["intermediate_periods"] = [200,100]
delete!(pf.config["method"]["options"]["ordering_error"], "ord_err_1") # remove from set
RPF.get_set_of_ordering_errors(pf) # Only ord_err_2 remains
pf.inputs[:ordering_error_functions]["ord_err_2"] = (
   (x,y) -> (sum((x[i] - y[i])^2 for i in eachindex(x)))
)

# Create and solve optimization problem
find_representative_periods(pf, reset=true)

# Save the values
u = pf.u; w = pf.w; v = pf.v; 
rep_periods = RPF.get_set_of_representative_periods(pf)

# Rebuild the periods finder and hot start it
config_file = normpath(joinpath(@__DIR__, "input_data", "default.yaml"))
pf = PeriodsFinder(config_file, populate_entries=true)
delete!(pf.config["method"], "clustering")
delete!(pf.config["method"]["optimization"], "time_series_error")
delete!(pf.config["method"]["options"], "ordering_error")
pf.config["method"]["optimization"]["duration_curve_error"]["type"] = "linear"
RPF.make_periods_finder_model!(pf)

# hot start values
for i in 1:length(u)
    set_start_value(pf.m.ext[:variables][:u][i], u[i])
    set_start_value(pf.m.ext[:variables][:w][i], w[i])
    for j in 1:length(rep_periods)
        jj = rep_periods[j]
        set_start_value(pf.m.ext[:variables][:v][i,jj], v[i,j])
    end
end

# optimize
optimize!(pf.m)
# TODO: function which assigns variables to pf.u, w, v, etc.+
create_plots(pf)
