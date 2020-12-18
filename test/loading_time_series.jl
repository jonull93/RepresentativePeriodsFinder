using RepresentativePeriodsFinder

config_file = joinpath(@__DIR__, "input_data", "default.yaml")
pf = PeriodsFinder(config_file)
pf.config["time_series"]["Missing Values"] = Dict(
    "source" => "missing_values.csv",
    "format"=> "yyyy-mm-dd HH:MM:SS",
    "value_column"=> "Missing values",
    "timestamp_column"=> "Timestamp",
)

ta = RepresentativePeriodsFinder.add_time_series!(pf, "Missing Values")

populate_entries!(pf)