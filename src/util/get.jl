# Sets
function get_set_of_time_series_names(pf::PeriodsFinder)
    return [k for (k,v) in pf.config["time_series"] if k != "default"]::Vector{String}
end

function get_set_of_time_series(pf::PeriodsFinder)
    return [k => v for (k,v) in pf.time_series if k != "default"]
end

function get_set_of_bins(pf::PeriodsFinder)
    return recursive_get(pf.config, 
        "method", "optimization", "duration_curve_error", "number_bins", 40
    )
end

function get_set_of_mandatory_periods(pf::PeriodsFinder)
    # TODO: fix this
    return Int64[]
end

function get_set_of_periods(pf::PeriodsFinder)
    return 1:get_number_of_periods(pf)
end

function get_set_of_representative_periods(pf::PeriodsFinder)
    if isdefined(pf, :u)
        return [p for p in u if p == 1]
    else
        return nothing
    end
end

function get_set_of_ordering_errors(pf::PeriodsFinder)
    opt_general = pf.config["method"]["options"]
    if haskey(opt_general, "ordering_error")
        return [k for (k,v) in opt_general["ordering_error"]]
    else
        return String[]
    end
end

# Simple parameters
function get_total_number_of_time_steps(pf::PeriodsFinder)
    opt = pf.config["method"]["options"]
    return opt["total_periods"] * opt["timesteps_per_period"]
end

function get_number_of_representative_periods(pf::PeriodsFinder)
    opt = pf.config["method"]["options"]
    return opt["representative_periods"]
end

function get_number_of_periods(pf::PeriodsFinder)
    opt = pf.config["method"]["options"]
    return opt["total_periods"]
end

function get_number_of_timesteps_per_period(pf::PeriodsFinder)
    opt = pf.config["method"]["options"]
    return opt["timesteps_per_period"]
end

function get_sampling_time(pf::PeriodsFinder)
    opt = pf.config["method"]["options"]
    return eval(Meta.parse(get(opt, "sampling_time", "Hour(1)")))
end

function get_error_term_weights(pf::PeriodsFinder)
    opt = pf.config["method"]["optimization"]
    opt_general = pf.config["method"]["options"]
    ord_errs = get_set_of_ordering_errors(pf)
    weights = Dict{String,Float64}()

    for err in ["time_series_error", "duration_curve_error"]
        weights[err] = recursive_get(opt, err, "weight", 0.0)
    end

    for err in ord_errs
        weights[err] = opt_general["ordering_error"][err]["weight"]
    end

    return weights
end

function get_time_series_weights(pf::PeriodsFinder)
    S = get_set_of_time_series(pf)
    return Dict(
        k => get(meta(ta), "weight", 0.0)
        for (k,ta) in S
    )
end

function get_ordering_error_function(pf::PeriodsFinder, err_name::String)
    opt_general = pf.config["method"]["options"]
    ord_err_opt = opt_general["ordering_error"]
    return errorfunc = 
        eval(Meta.parse(ord_err_opt[err_name]["function"], raise=true))
end

# Parameters
function get_normalised_time_series_values(pf::PeriodsFinder)
    if isdefined(pf, :x) == false || isempty(pf.x) == true
        for (ts_name, ta) in get_set_of_time_series(pf)
            ts_name = meta(ta)["name"]
            start = meta(ta)["start"]

            sampling_time = get_sampling_time(pf)
            ntt = get_total_number_of_time_steps(pf)
            timestamps = range(start, length=ntt, step=sampling_time)
            vec = values(ta[timestamps])
            @assert length(vec) == ntt "Less than $ntt values lie between $(timestamps[1]) and $(timestamps[end]) in $ts_name"

            vec = normalize_values(vec)

            npt = get_number_of_periods(pf)
            ntpp = get_number_of_timesteps_per_period(pf)
            pf.x[ts_name] = permutedims(reshape(vec, ntpp, npt), (2,1))
        end
    end
    return pf.x
end

function get_duration_curve_parameter(pf::PeriodsFinder)

    bins = get_set_of_bins(pf)
    periods = get_set_of_periods(pf)
    A = Dict(ts_name => Array{Float64,2}(undef, length(bins), length(periods)))
    for (ts_name, ta) in get_set_of_time_series(pf)
        cum_bin_end!(A, bins, periods, ta)
    end

    # ##################################################################################
    # # Calculate bins and matrices
    # ##################################################################################
    # self.matrix_full = reshape(self.data, round(Int, length(self.data) / length(dft.periods)), length(dft.periods))'

    # normalise_time_series!(self)

    # value_width = (maximum(self.data_norm) - minimum(self.data_norm)) / length(dft.bins)

    # mini = minimum(self.data_norm)
    # bins = [mini; [mini + i * value_width for i in range(1,stop=length(dft.bins))]]
    # bins[end] += 0.01

    # self.matrix_bins = Array{Int}(undef, length(dft.periods), length(bins)-1)

    # for p in 1:size(self.matrix_full_norm)[1]
    #     self.matrix_bins[p,:] = fit(Histogram, self.matrix_full_norm[p,:], bins, closed=:left).weights'
    # end

    # self.matrix_bins_cumsum = cumsum(sum(self.matrix_bins,dims=1)/length(self.data_norm),dims=2)[:] # -> L

    # self.matrix_bins = self.matrix_bins / round(Int,length(self.data_norm)/length(dft.periods))

    # self.matrix_bins_cumsum_day = cumsum(self.matrix_bins, dims=2) # -> A
end