using RepresentativePeriodsFinder

config_file = joinpath(@__DIR__, "input_data", "default.yaml")
pf = PeriodsFinder(config_file)

# Add and process load (easy)
ta = RepresentativePeriodsFinder.read_time_series(pf, "Load")
ta = RepresentativePeriodsFinder.interpolate_missing_values!(pf, ta)
ta = RepresentativePeriodsFinder.resample!(pf, ta)

ta = RepresentativePeriodsFinder.read_time_series(pf, "Residual Load")

# ta = RepresentativePeriodsFinder.adjust_length!(pf, ta)

pf.config["time_series"]["Missing Values"] = Dict(
    "source" => "missing_values.csv",
    "format"=> "yyyy-mm-dd HH:MM:SS",
    "value_column"=> "Missing values",
    "timestamp_column"=> "Timestamp",
)
ta = RepresentativePeriodsFinder.read_time_series(pf, "Missing Values")

pf.config["time_series"]["Linear Increase to Resample"] = Dict(
    "source" => "linear_increasing_24_hours.csv",
    "value_column"=> "Value",
    "timestamp_column"=> "Timestamp",
    "start" => "",
    "delim" => ','
)
ta = RepresentativePeriodsFinder.read_time_series(pf, "Linear Increase to Resample")
ta = RepresentativePeriodsFinder.interpolate_missing_values!(pf, ta)
ta = RepresentativePeriodsFinder.resample!(pf, ta)

populate_entries!(pf)


# Create a linearly increasing value with semi random time stamps
# timestamps = DateTime(1970,1,1) .+ [
#     Hour(0),
#     Hour(1),
#     Hour(3) + Minute(15),
#     Hour(7) + Minute(30),
#     Hour(20) + Minute(10),
#     Hour(20) + Minute(15),
#     Hour(20) + Minute(40),
#     Hour(24) + Minute(30)
# ]
# vals = Dates.value.(timestamps)

# y = (vals .- vals[1])./(vals[end] - vals[1])
# df = DataFrame(
#     Timestamp = timestamps,
#     Value = y
# )
# CSV.write(joinpath(@__DIR__, "input_data", "linear_increasing_24_hours.csv"), df)