# Loading time series
In the absence of additional documentation, refer to the examples in the [`test`](https://gitlab.kuleuven.be/UCM/representativedaysfinder.jl/-/tree/dev/test) directory of this repository to understand the options available for loading in time series.

For the purposes of illustration, consider [`default.yaml`](https://gitlab.kuleuven.be/UCM/representativedaysfinder.jl/-/blob/dev/test/input_data/default.yaml) which is used for testing this package.

```yaml
time_series:
  default: # These values will be used unless specified in the entries below
    source: "time_series.csv"
    csv_options: # These are passed to CSV.read
        delim: ";"
    timestamp: "Timestamp" # optional
    weight: 1 # weight of time series in objective / clustering
    silence_warnings: false # If true, won't throw error msgs when reading csv
    interpolation_type: "linear" # or "constant" or don't specify
    sampling_time: "Hour(1)" # Only necessary if a timestamp is not given
  
  Residual Load:
    source: "residual_load.csv"
    csv_options: 
      dateformat: "yyyy-mm-dd HH:MM:SS"
      # delim: ";"
      missingstring: "#N/A"
      # missingstrings: ["NA", "N/A", "na"] # Careful with the plural!
      silencewarnings: true
    value_column: "Residual Load [MW]"
    timestamp_column: "Timestamp"
    
    start: "2017-01-01 00:00:00"
    resample: true # or false or don't specify

  Load:
    value_column: "Load"
    start: 1 # Can also specify an index in the CSV file (not including header)

  LFW:
    value_column: "LFW"
    delim: ';'
    weight: 1

  LFS:
    value_column: "LFS"
    weight: 1
```

To load the above time series (without selecting representative periods) run [`pf = PeriodsFinder(config_file)`](@ref) where `config_file` is the path to your chosen configuration file.

The `default` entry specifies options to be used if these are not specified for the time series entries. So `Load` is read in from `time_series.csv` since there is no `source` option specified for it, while `Residual Load` is read from `residual_load.csv`.