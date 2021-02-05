# Outputs

## Saving and inspecting results
If the option `save_results` or `create_plots` in `results` of the configuration file is true, then running [`find_representative_periods`](@ref) will save the selection results to `.csv` and create duration curve plots in the specified directory. 

Otherwise the following functions can be called:
```julia
dir = "/home/user/Desktop/selecting_periods/results"
save(pf, dir)
create_plots(pf, dir)
```

If no directory is specified then the directory specified in `results: result_dir` will be used.

## Loading results
Results saved using [`save(pf, dir)`](@ref) can also be loaded using [`load(pf, dir)`](@ref).