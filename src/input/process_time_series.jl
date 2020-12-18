function interpolate_missing_values!(pf::PeriodsFinder, ts_name::String)
    ta = pf.time_series[ts_name]
    ts_dict = pf.config["time_series"]
    interpolation_type = config_get(
        ts_dict[ts_name], ts_dict["default"], "interpolation_type", "linear"
    )
    itp = make_interpolator(ta, interpolation_type)
    A = values(ta)
    
    n_missing = sum(ismissing.(A))
    @debug "Interpolating for $n_missing values"

    for i in 1:length(ta)
        val = values(ta[i])[1]
        x = time_to_grid_val(timestamp(ta[1])[1], timestamp(ta[i])[1])
        if ismissing(val)
            A[i,:] .= itp(x)
        end
    end
end

function make_interpolator(ta::TimeArray, interpolation_type::String)
    y = eltype(values(ta))[]
    grid = Int64[]
    for i in 1:length(ta)
        val = values(ta[i])[1]
        x = time_to_grid_val(timestamp(ta[1])[1], timestamp(ta[i])[1])
        if ismissing(val) == false
            push!(y, val)
            push!(grid, x)
        end
    end

    if interpolation_type == "linear"
        itp = interpolate((grid,), y, Gridded(Linear()))
        itp = extrapolate(itp, Flat())
        return itp
    elseif interpolation_type == "constant"
        itp = interpolate((grid,), y, Gridded(Constant()))
        itp = extrapolate(itp, Flat())
        return itp
    end

end

function time_to_grid_val(t_init::DateTime, t::DateTime)
    return (t - t_init).value
end