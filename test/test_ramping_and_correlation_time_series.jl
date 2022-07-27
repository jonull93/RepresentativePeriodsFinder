# # Ramping and correlation time series

# In the original paper by [Poncelet et al.](https://www.mech.kuleuven.be/en/tme/research/energy_environment/Pdf/wp-2015-10b.pdf) there was the possibility of using a correlation time series, e.g. for electrical load and solar power (which are highly correlated over the course of a day). Equally, it may be interesting to select representative periods based on a ramping or mileage time series (e.g. to select days with a high rate of change in the load for example).

# Eventually, this could be added as a core feature within `RepresentativePeriodsFinder.jl`, and indeed there are still two open issues regarding this (see [here](https://gitlab.kuleuven.be/UCM/representativeperiodsfinder.jl/-/issues/16) and [here](https://gitlab.kuleuven.be/UCM/representativeperiodsfinder.jl/-/issues/17)). This is not a priority however and will not be done unless someone (that could be you!) opens a pull request. 

# Anyway, enough faffing about. First, the ramping time series (since it's easier):