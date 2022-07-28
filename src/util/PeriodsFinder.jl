"""
    PeriodsFinder(config_file::String; populate_entries::Bool=false)

Type used for finding representative periods in a time series (see `find_representative_periods`). If `populate_entries=true` then time series are loaded into `pf`.

# Entries

* `config_file`: Path of configuration file used to create the `PeriodsFinder`.
* `config`: Dictionary of input parameters.
* `time_series`: Dictionary of time series.
* `x`: Normalise time series in matrix format (rows = timesteps, columns = periods).
* `inputs`: Used internally to save variables which are computationally intensive to calculate.
* `m`: JuMP model used for representative period selection.
* `u`: Selection variable, a vector whose entries are 1 if a period is selected and 0 otherwise.
* `w`: Weighting variable, a vector whose entries greater than 0 if a period is selected and 0 otherwise.
* `v`: Ordering variable, a matrix where the sum over the columns is equal to `w` and the diagonal is equal to `u`.  
"""
mutable struct PeriodsFinder
    ###########################################################################
    # Configuration
    ###########################################################################
    config_file::AbstractString
    config::Dict{Any,Any}

    ###########################################################################
    # Inputs
    ###########################################################################
    time_series::Dict{String,FloatTimeArray}
    x::Dict{String,Array{Float64,2}} # Normalised time series in matrix format
    inputs::Dict{Symbol,Any} # Any saved inputs (sets, parameters, etc)
    # TODO: Change this to Dict{Union{Symbol,String},Any} - much nicer.

    ###########################################################################
    # optimization model
    ###########################################################################
    m::JuMP.Model

    ###########################################################################
    # Results
    ###########################################################################
    u::Array{Bool,1}    # Selection variable
    w::Array{Float64,1} # Weight variable (can also be integer)
    v::Array{T,2} where T <: Union{Bool,Float64} # Ordering variable (not always defined, can be float)

    function PeriodsFinder()
        pf = new()
        pf.config = Dict{Any,Any}()
        pf.time_series = Dict{String,TimeArray}()
        return pf
    end

    function PeriodsFinder(config_file::String;     
            populate_entries::Bool=false
        )
        pf = new()
        pf.config_file = config_file
        @assert isfile(config_file)
        pf.config = YAML.load(open(config_file))
        !(haskey(pf.config, "base_dir")) && (pf.config["base_dir"] = dirname(config_file))
        pf.time_series = Dict{String,TimeArray}()

        if populate_entries == true
            pf = populate_entries!(pf)
        end

        return pf
    end

    # Pretty prints of long term planning model
    function Base.print(io::IO, pf::PeriodsFinder)
        println(io, "Representative Periods Finder")
        isdefined(pf, :config_file) && println(io, "Configuration file: \n\t$(pf.config_file)")
    end

    function Base.show(io::IO, pf::PeriodsFinder)
        print(io, pf)
    end
end

"""
    populate_entries!(pf::PeriodsFinder; reset_inputs::Bool=true)

Re-reads time series from the `.csv` files. If `reset_inputs=true`, also calls `reset_inputs!(pf)`.
"""
function populate_entries!(pf::PeriodsFinder)
    S = get_set_of_time_series_names(pf)
    
    for ts_name in S
        @info "Adding $ts_name..."
        ta = read_time_series(pf, ts_name)
        if haskey(meta(ta), "interpolation_type")
            @info "Interpolating missing values..."
            interpolate_missing_values!(pf, ta)
        end
        if haskey(meta(ta), "resample") && meta(ta)["resample"] == true
            @info "Resampling..."
            resample!(pf, ta)
        end
        pf.time_series[ts_name] = TimeArray(
            DateTime.(timestamp(ta)), Float64.(values(ta)), 
            colnames(ta), meta(ta)
        )
    end

    reset_inputs!(pf)

    return pf
end

"""
    reset_inputs!(pf::PeriodsFinder)

Resets the inputs of `pf`, specifically `pf.x` and `pf.inputs`. Should be called after making changes to `pf.config` to ensure that these are applied when selecting representative periods.

# Example

```jldoctest
pf = PeriodsFinder(RepresentativePeriodsFinder.datadir("only_load.yaml"); populate_entries=true);
println("Length of L: $(length(RepresentativePeriodsFinder.get_discretised_duration_curve(pf)["Load"]))")

pf.config["method"]["optimization"]["duration_curve_error"]["number_bins"] = 80;
println("Length of L: $(length(RepresentativePeriodsFinder.get_discretised_duration_curve(pf)["Load"]))")

reset_inputs!(pf);
println("Length of L: $(length(RepresentativePeriodsFinder.get_discretised_duration_curve(pf)["Load"]))")

# output

[ Info: Adding Load...
[ Info: Interpolating missing values...
Length of L: 40
Length of L: 40
Length of L: 80
```
"""
function reset_inputs!(pf::PeriodsFinder)
    pf.x = Dict{String,Array{Float64,2}}()

    # Instantiate inputs, load error functions
    ord_err_names = get_set_of_ordering_errors(pf)
    pf.inputs = Dict(
        :ordering_error_functions => Dict{String,Function}(
            ord_err => get_ordering_error_function(pf, ord_err)
            for ord_err in ord_err_names
            # Do this now to avoid "new world" issues
        )
    )

    return pf
end