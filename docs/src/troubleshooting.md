# Troubleshooting

This package is tested using Julia 1.2 and above. If it doesn't work for you, this is in all likelihood due to your own environment! To fix this, start from a blank environment and add this package:

```julia
using Pkg
Pkg.activate(".") # Create a new environment in this directory, which should not contain a Project.toml file
Pkg.add(url="https://gitlab.kuleuven.be/UCM/representativedaysfinder.jl")
```

See also the troubleshooting tips below.

## GR

If issues with the GR-engine (i.e. plotting) occur then build the `GR` package again:

```julia
]build GR
```

## CSV

Many versions of the `CSV` package are incompatible with each other, so errors which mention `CSV` are likely due to some other package restricting later versions of `CSV` from being used. This is a you problem unfortunately, so you'll have to fix it yourself (by e.g. removing the troublesome package).

## DataFrames

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

## "Key not found"

If you get an error that mentions this, forcing the `PeriodsFinder` type to load your time series with `populate_entries=true` may help:

```julia
config_file = <path_to_config_file>
pf = PeriodsFinder(config_file; populate_entries=true) # This keyword argument is the "key" (ha)
find_representative_periods(pf; optimizer=Cbc.Optimizer)
```

## Optimal ordering of representative periods is slow

It may be that setting the `integral_weights` option to `true` speeds up the optimisation massively (from several hours or days to several seconds).

## Other errors and fixes

* Make sure your time series are correct, i.e. without `"NaN"`, `"N/A"` or blank entries.
* The time series interpolation feature could fix `"NaN"` entries, but it is quite experimental so best to turn it off.
* If you've specified a timestamp column for a time series, not specifying it could fix possible bugs.
* The most bulletproof workflow is to specify all your options in the configuration file and then call [`find_representative_periods`](@ref). Modifying `pf.config` and then calling [`find_representative_periods`](@ref)] may have unintented consequences. If you have to do it, remember to call [`populate_entries!`](@ref) or [`reset_inputs!`](@ref) beforehand.
* This package is not compatible with 32 bit machines or Julia versions. Then again, why on earth are would you do that?
* For "age of world" errors, see [Selecting representative periods](@ref).
