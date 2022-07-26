using RepresentativePeriodsFinder, JuMP, Cbc
RPF = RepresentativePeriodsFinder

# Create PeriodsFinder
config_file = RPF.datadir("default.yaml")
pf = PeriodsFinder(config_file, populate_entries=true)

# Delete entries for time series and duration curve error
delete!(pf.config["method"]["optimization"], "time_series_error")
delete!(pf.config["method"]["options"], "ordering_error")
delete!(pf.config["method"], "clustering")

# Choose number of rep hours
n_rep = 5
pf.config["method"]["options"]["representative_periods"] = n_rep

# Make life a bit easier for the optimiser, reduce the number of hours
pf.config["method"]["options"]["total_periods"] = 20

# Change period length to 1 hour per period
pf.config["method"]["options"]["time_steps_per_period"] = 1

# Setup optimizer
opt = optimizer_with_attributes(Cbc.Optimizer, "seconds" => 60)

# Use absolute error (this is the only possibility with Cbc)
pf.config["method"]["optimization"]["duration_curve_error"]["type"] = "absolute"
find_representative_periods(pf, optimizer=opt, reset=true)

# Assert that the ordering variable is not made
@testset "Optimisation based hour selection" begin
    @test typeof(pf.m.ext[:variables][:v]) <: RPF.SingleValuedContainer
    @test haskey(pf.m.ext[:constraints], :restrict_ordering) == false
    @test haskey(pf.m.ext[:constraints], :linear_combination_representative_periods) == false
    @test length(RPF.get_set_of_representative_periods(pf)) == n_rep
end