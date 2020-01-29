##################################################################################
# Author:   Hanspeter Höschle
# Date:     15/06/2017
##################################################################################

##################################################################################
# Create basic TimeSeries from config-file
##################################################################################
# import RepresentativeDaysFinders: TimeSeries
# using RepresentativeDaysFinders: DaysFinderTool, calculate_matrix_bins!
function TimeSeries(dft::DaysFinderTool, config::Dict{Any,Any})
    self = TimeSeries()
    self.config = config
    self.name = config["name"]
    self.time_series_type = :basic
    self.weight = self.config["weight"]

    ##################################################################################
    # Load data
    ##################################################################################
    if "csv" in keys(config["source"])
        col = config["source"]["column"]
        csv_file = normpath(joinpath(dft.config["basedir"], config["source"]["csv"]))
        df = CSV.read(csv_file, delim = config["source"]["delimiter"][1])
        self.data = df[!,Symbol(col)]
    end

    calculate_matrix_bins!(self, dft)

    return self
end
##################################################################################
# Create diff TimeSeries from existing TimeSeries
##################################################################################
function TimeSeries(dft::DaysFinderTool, basic::TimeSeries, method::Symbol=:all_1step)
    self = TimeSeries()
    self.time_series_type = :diff
    self.weight = 1

    if method == :all_1step
        self.name = "diff1_$(basic.name)"
        self.data = [0.; diff(basic.data)]
    else
        error()
    end

    calculate_matrix_bins!(self, dft)

    return self
end

##################################################################################
# Create correlation TimeSeries from two TimeSeries
##################################################################################
function TimeSeries(dft::DaysFinderTool, ts1::TimeSeries, ts2::TimeSeries)
    self = TimeSeries()
    self.time_series_type = :diff
    self.weight = 1
    self.name = "corr_$(ts1.name)_$(ts2.name)"

    # Normalization
    data1 = ts1.data / sum(ts1.data) * 10000
    data2 = ts2.data / sum(ts2.data) * 10000

    # Time series based on numerator correlation definition
    self.data = (data1 .- mean(data1)) .* (data2 .- mean(data2))

    calculate_matrix_bins!(self, dft)

    return self
end

function calculate_matrix_bins!(self::TimeSeries, dft::DaysFinderTool)
    ##################################################################################
    # Calculate bins and matrices
    ##################################################################################
    self.matrix_full = reshape(self.data, round(Int,length(self.data) / length(dft.periods)), length(dft.periods))'

    value_width = (maximum(self.data) - minimum(self.data)) / length(dft.bins)

    mini = minimum(self.data)
    bins = [mini; [mini + i * value_width for i in range(1,stop=length(dft.bins))]]
    bins[end] += 0.01

    self.matrix_bins = Array{Int}(undef, length(dft.periods), length(bins)-1)

    for p in 1:size(self.matrix_full)[1]
        self.matrix_bins[p,:] = fit(Histogram, self.matrix_full[p,:], bins, closed=:left).weights'
    end

    self.matrix_bins_cumsum = cumsum(sum(self.matrix_bins,dims=1)/length(self.data)*100.,dims=2)[:] # -> L

    self.matrix_bins = 100. * self.matrix_bins / round(Int,length(self.data)/length(dft.periods))

    self.matrix_bins_cumsum_day = cumsum(self.matrix_bins, dims=2) # -> A
end

function weight!(w, ts::TimeSeries)
    w[ts.name] = ts.weight
end

function area_total!(at, ts::TimeSeries)
    at[ts.name] = sum(abs.(ts.data))
end

function area_total_days!(atd, periods, ts::TimeSeries)
    for (idx, p) in enumerate(periods)
        atd[ts.name, p] = sum(abs.(ts.matrix_full[idx, :]))
    end
end

function cum_bin_end!(a, periods, bins, ts::TimeSeries)
    for (idx_p, p) in enumerate(periods)
        for (idx_b, b) in enumerate(bins)
            a[ts.name, p, b] = ts.matrix_bins_cumsum_day[idx_p, idx_b]
        end
    end
end

function cum_bin_total!(l, bins, ts::TimeSeries)
    for (idx_b, b) in enumerate(bins)
        l[ts.name, b] = ts.matrix_bins_cumsum[idx_b]
    end
end
