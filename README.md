# RepresentativePeriodsFinder

This is a Julia package to select representative periods from time series data, typically to then be used in capacity expansion planning models. For more information see the [documentation]().

## Installation

* `RepresentativePeriodsFinder` can be added via the Julia package manager (type `]`): `pkg> add https://gitlab.kuleuven.be/UCM/representativedaysfinder.jl`

## Basic useage

```julia
julia> using RepresentativePeriodsFinder, JuMP, Cbc
julia> config_file = <path_to_config_file>
julia> find_representative_periods(config_file, optimizer=Cbc.Optimizer)
```
For more information see the [documentation]().

## Trouble shooting
If issues with GR-engine occur build the `GR` package:

```julia
pkg> build GR
```

## Reporting issues and feature requests

If you have any issues or feature requests, report them on GitLab and mention `@steffen.kaminski` of `@u0128861`. Most of the GitLab issues already there are just waiting for someone to motivate their implementation.

## Alternative packages

* [TimeSeriesClustering](https://holgerteichgraeber.github.io/TimeSeriesClustering.jl/stable/quickstart/)

## Developers

```julia
julia> ENV["JULIA_PKG_DEVDIR"] = "<path_to_dir_of_choice>"
pkg> dev "https://gitlab.kuleuven.be/UCM/representativedaysfinder.jl"
```


## License

...
