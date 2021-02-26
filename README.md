# RepresentativePeriodsFinder

This is a Julia package to select representative periods from time series data, typically to then be used in capacity expansion planning models. For more information see the [documentation](https://ucm.pages.gitlab\.kuleuven\.be/representativeperiodsfinder.jl/).

## Installation

* `RepresentativePeriodsFinder` can be added via the Julia package manager (type `]`): `pkg> add https://gitlab.kuleuven.be/UCM/representativedaysfinder.jl`

## Basic useage

```julia
using RepresentativePeriodsFinder, Cbc
config_file = <path_to_config_file>
find_representative_periods(config_file, optimizer=Cbc.Optimizer)
```
For more information see the [documentation](https://ucm.pages.gitlab\.kuleuven\.be/representativeperiodsfinder.jl/).

## Trouble shooting

### GR
If issues with GR-engine occur build the `GR` package:

```julia
pkg> build GR
```

### CSV
For some reason the `CSV` package can be troublesome, so make sure that this package is updated in your main environment by doing:
```julia
pkg> add CSV
pkg> up CSV
```

## Reporting issues and feature requests

If you have any issues or feature requests, report them on GitLab and mention `@steffen.kaminski` or `@u0128861`. Most of the GitLab issues already there are just waiting for someone to motivate their implementation.

## Alternative packages

* [TimeSeriesClustering.jl](https://holgerteichgraeber.github.io/TimeSeriesClustering.jl/stable/quickstart/)

## Developers

```julia
julia> ENV["JULIA_PKG_DEVDIR"] = "<path_to_dir_of_choice>"
pkg> dev "https://gitlab.kuleuven.be/UCM/representativedaysfinder.jl"
```

## Acknowledgements
This package is largely based on the [paper](https://www.mech.kuleuven.be/en/tme/research/energy_environment/Pdf/wp-2015-10b.pdf) and work of Kris Poncelet and Hanspeter Höschle where optimisation instead of clustering techniques were used to select representative days. This method was initially implemented in GAMS. An first Julia version was later implemented by Hanspeter Höschle while he was working at [VITO](https://vito.be/en). This was further developed by Sebastian Gonzato for a [paper on long term storage](https://www.mech.kuleuven.be/en/tme/research/energy-systems-integration-modeling/pdf-publications/wp-esim2021-1).

## License

This project is licensed under the GNU General Public License (GPL) v3.0 or later - see the `LICENSE.md` file for details.
