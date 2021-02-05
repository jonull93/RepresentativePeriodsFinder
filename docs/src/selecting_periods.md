# Selecting representative periods
In the absence of additional documentation, refer to the examples in the [`test`](https://gitlab.kuleuven.be/UCM/representativedaysfinder.jl/-/tree/dev/test) directory of this repository to understand options availabile when selecting representative periods.

To select periods, the general workflow is:
```julia
using RepresentativePeriodsFinder, JuMP, Cbc
config_file = "/home/user/Desktop/selecting_periods/config_file.yaml"
pf = PeriodsFinder(config_file)
optimizer = optimizer_with_attributes(Cbc.Optimizer, "seconds" => 300) # seconds specifies time out limit
pf = find_representative_periods(pf, optimizer=optimizer)
```

Alternatively, you can just run:
```julia
config_file = "/home/user/Desktop/selecting_periods/config_file.yaml"
optimizer = optimizer_with_attributes(Cbc.Optimizer, "seconds" => 300)
pf = find_representative_periods(config_file, optimizer=optimizer)
```

You do not need to specify an optimizer if you're using a clustering based method.

Running [`find_representative_periods(pf)`](@ref) will produce three outputs:
* `pf.u` - a boolean vector indicating which periods have been selected to be representative.
* `pf.w` - a vector of the weights attached to periods (0 if thse are not representative).
* `pf.v` - a matrix mapping representative periods (columns) to non-representative periods (rows). Depending on the selection method chosen, this may not be defined.

## Common options

## Optimisation methods
There is one general optimisation problem which is solved with differences coming from the objective terms and whether periods are selected or also ordered throughout the year.

## Clustering
You must specify your own dissimilarity function(s) in `options: ordering_error`.