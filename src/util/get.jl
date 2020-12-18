function get_set_of_time_series(pf::PeriodsFinder)
    return [k for (k,v) in pf.config["time_series"] if k != "default"]
end