function get_set_of_time_series(pf::PeriodsFinder)
    return [k for (k,v) in pf.config["time_series"] if k != "default"]
end

function get_total_number_of_time_steps(pf::PeriodsFinder)
    return pf.config["total_periods"] * pf.config["timesteps_per_period"]
end