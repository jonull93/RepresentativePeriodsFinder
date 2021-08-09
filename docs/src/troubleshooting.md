# Troubleshooting

This package is tested using Julia 1.2, 1.3, 1.4 and 1.5. If it doesn't work for you, this is in all likelihood due to your own environment! Some troubleshooting tips to fix this are suggested below.

## GR

If issues with the GR-engine (i.e. plotting) occur build the `GR` package:

```julia
pkg> build GR
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

* Make sure your time series are correct, i.e. without "NaN", "N/A" or blank entries.
* The time series interpolation feature could fix "NaN" entries, but it is quite experimental so best to turn it off.
* If you've specified a timestamp column for a time series, not specifying it could fix possible bugs.
