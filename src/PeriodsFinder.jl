###############################################################################
# Author:   Hanspeter Höschle
# Date:     15/06/2017
###############################################################################
mutable struct PeriodsFinder
    ###########################################################################
    # Configuration
    ###########################################################################
    config_file::AbstractString
    config::Dict{Any,Any}

    ###########################################################################
    # Inputs
    ###########################################################################
    time_series::Dict{String,TimeArray}
    x::Dict{String,Array{Float64,2}} # Normalised time series in matrix format
    inputs::Dict{Any,Any} # Any other inputs (sets, parameters, etc)

    ###########################################################################
    # Sets
    ###########################################################################
    # bins::Array                     # Set of bins
    # curves::Array{String,1}         # Set of duration curves
    # periods::Array{Int64,1}         # Set of potential representative periods
    # rep_periods::Array{Int64,1}     # Set of chosen representative periods - list of integers
    # timesteps::UnitRange            # Set of timesteps per potential period

    ###########################################################################
    # Parameters
    ###########################################################################
    # N_representative_periods::Int   # Number of representative periods to select
    # N_total_periods::Int            # Total number of repetitions required 
                                      # to scale up the duration of a single 
                                      # representative period to one year

    # WEIGHT_DC::Dict                         # Relative importance of approximating duration curve c
    A::Dict                                 # Share of the time of day d during which the lowest value of the range corresponding to bin b of duration curve c is exceeded
    L::Dict                   # Share of the time during which the values of a time series with corresponding duration curve c exceed the lowest value of the range corresponding to bin b
    # AREA_TOTAL::Dict                        # Total area under the DC of time series c
    # AREA_TOTAL_DAY::Dict      # Total area under time series with DC c per day
    # AREA_ERROR_TOLERANCE::Dict              # tolerance parameter for area constraint

    ###########################################################################
    # Optimisation model
    ###########################################################################
    m::JuMP.Model

    ###########################################################################
    # Results
    ###########################################################################
    u::Array{Int64,1} # Selection variable
    w::Array{Float64,1} # Weight variable
    v::Array{Float64,2} # Ordering variable (not always defined)

    function PeriodsFinder(config_file::String; populate_entries::Bool=false)
        pf = new()
        pf.config_file = config_file
        pf.config = YAML.load(open(config_file))
        pf.config["base_dir"] = dirname(config_file)
        pf.time_series = Dict{String,TimeArray}()

        if populate_entries == true
            pf = populate_days_finder!(pf)
        end

        return pf
    end

    # Pretty prints of long term planning model
    function Base.print(io::IO, dft::PeriodsFinder)
        println(io, "Representative Periods Finder")
        println(io, "Configuration file: \n\t$(dft.config_file)")
    end

    function Base.show(io::IO, dft::PeriodsFinder)
        print(io, dft)
    end
end

function populate_days_finder!(self::PeriodsFinder)
    pad = 3

    ##################################################################################
    # Initialize the TimeSeries
    ##################################################################################
    self.time_series = Dict()
    self.curves = Array{String,1}()
    for ts_config in self.config["time_series"]
        println("-"^100)

        ts = TimeSeries(self, ts_config)
        addTimeSeries!(self, ts)

        println("$(ts.name) added")

        ##################################################################################
        # Dynamics profiles
        ##################################################################################
        if self.config["dynamics_method"] == "all_1step"
            println("-"^100)

            ts_diff = TimeSeries(self, ts, :all_1step)
            addTimeSeries!(self, ts_diff)

            println("$(ts_diff.name) added")
        end
    end

    ##################################################################################
    # Parameters
    ##################################################################################
    self.WEIGHT_DC                  = Dict()
    self.A                          = Dict()
    self.L                          = Dict()

    self.AREA_TOTAL                 = Dict()
    self.AREA_TOTAL_DAY             = Dict()

    self.N_representative_periods   = try_get_val(self.config, "number_days", 8)
    self.N_total_periods            = try_get_val(self.config, "number_days_total", 365)

    ##################################################################################
    # Sets
    ##################################################################################
    self.bins =     ["b" * string(b, pad=pad) for b in range(1, stop=self.config["number_bins"])]

    lenP = self.N_total_periods
    self.periods = 1:lenP

    # TODO: this should just be timeseries length / N_periods
    # But this makes an assumption on the length of the timeseries, so eh...
    lenT = try_get_val(self.config, "timesteps_per_period", 24)
    self.timesteps = 1:lenT

    ##################################################################################
    # Correlation
    ##################################################################################
    if self.config["correlation_method"] == "all"
        ts_basic = [ts.name for ts in values(self.time_series) if ts.time_series_type == :basic]
        for combi in combinations(ts_basic, 2)
            println("-"^100)
            ts1 = self.time_series[combi[1]]
            ts2 = self.time_series[combi[2]]
            ts = TimeSeries(self, ts1, ts2)
            addTimeSeries!(self, ts)

            println("$(ts.name) added")
        end
    end

    # misc dictionary
    self.misc = Dict()

    return self
end