# RepresentativePeriodsFinder

This is a Julia package to select representative periods from time series data, typically to then be used in capacity expansion planning models. For more information see the [documentation](https://ucm.pages.gitlab\.kuleuven\.be/representativeperiodsfinder.jl/).

## Installation

* Download and install Julia from [here](https://julialang.org/downloads/).
* Assuming you've added Julia to your path, you can open up a Julia REPL by typing `julia` in a terminal.
* `RepresentativePeriodsFinder` can be added via the Julia package manager (type `]`): `pkg> add https://gitlab.kuleuven.be/UCM/representativedaysfinder.jl`

## Basic useage

```julia
using RepresentativePeriodsFinder, Cbc
config_file = <path_to_config_file>
find_representative_periods(config_file, optimizer=Cbc.Optimizer)
```
For more information see the [documentation](https://ucm.pages.gitlab\.kuleuven\.be/representativeperiodsfinder.jl/).

## Troubleshooting

This package is tested using Julia 1.2, 1.3, 1.4 and 1.5. If it doesn't work for you, this is in all likelihood due to your own environment! Some troubleshooting tips to fix this are suggested below.

### GR
If issues with GR-engine occur build the `GR` package:

```julia
pkg> build GR
```

### CSV
Many versions of the `CSV` package are incompatible with each other, so errors which mention `CSV` are likely due to some other package restricting later versions of `CSV` from being used. This is a you problem unfortunately, so you'll have to fix it yourself (by e.g. removing the troublesome package).

### DataFrames
`DataFrames.jl` can also be quite a troublesome package, since indexing changed between versions. If you get an error about a column not being in a `DataFrame` (i.e. your `.csv` file) when you know it should be, make sure your version of `DataFrames` is at least `0.20.0`:
```julia
]add DataFrames
]status # Check what version number is written next to the entry with DataFrames
```
If it's not, try to update:
```julia
]up DataFrames
```
If this doesn't work, then it's likely some other package is preventing you from updating.

### "Key not found"
If you get an error that mentions this, forcing the `PeriodsFinder` type to load your time series may help:
```julia
config_file = <path_to_config_file>
pf = PeriodsFinder(config_file; populate_entries=true) # This keyword argument is the "key" (ha)
find_representative_periods(pf; optimizer=Cbc.Optimizer)
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
