
config_file = normpath(joinpath(pwd(), "..", "scenarios", "2017_interconnected_cm_2015.yaml"))
config = YAML.load(open(config_file))
dft = RepresentativeDaysFinder.DaysFinderTool(config_file)

ts_config = config["time_series"][1]

ts = RepresentativeDaysFinder.TimeSeries(dft, ts_config)
periods = RepresentativeDaysFinder.get_mandatory_periods(ts, dft)