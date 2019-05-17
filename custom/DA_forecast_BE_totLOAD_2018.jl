using RepresentativeDaysFinders
using JuMP
using PyCall
# using GLPK
using Gurobi
using CSV
using DataFrames
using Dates
using TimeZones

##################################################################################
# collect data
##################################################################################
#done externally

##################################################################################
# preprocess data if needed
##################################################################################
#interpolation done here
# df_LF_wind = ENTSOEcsv2dataframe(joinpath(@__DIR__, "windforecast_ELIA.csv"), 12, :LF_wind; delim = ';')
df_LF_wind = CSV2DataFrame(joinpath(@__DIR__, "windforecast_ELIA.csv"); delim = ';')
df_LF_wind_filled = RepresentativeDaysFinders.check_temporal_consistency(df_LF_wind, 0.25, "y-m-d HH:MM:SSz", :Column1)
df_LF_wind_interpolated = interpolatedataframe(df_LF_wind_filled,:DateFormated,:LoadFactor, 1.)

# itp = RepresentativeDaysFinders.makeinterpolator(df_LF_wind_filled, :LoadFactor)
df_LF_solar = ENTSOEcsv2dataframe(joinpath(@__DIR__, "solarforecast_ELIA.csv"), 2, :DateFormated; delim = ';')
# df_LF_solar_interpolated = interpolatedataframe(df_LF_solar, :LF_solar, 0.25,1.)
# csv_wind_raw = CSV.read(joinpath(@__DIR__, "windforecast_ELIA.csv"), delim = ";")
# csv_solar_raw = CSV.read(joinpath(@__DIR__, "solarforecast_ELIA.csv"), delim = ";")
# df_test = join(csv_wind_raw, csv_solar_raw, on = [:Column1, :Column1], makeunique=true )





test2 = CSV2DataFrame(joinpath(@__DIR__, "windforecast_ELIA.csv"), 1, :bla; delim = ';')

##################################################################################
# Specify location of config-file
##################################################################################
using Cbc
config_file = joinpath(@__DIR__,"DA_LOAD_forcast_BE_2018.yaml")
dft = findRepresentativeDays(config_file, with_optimizer(Cbc.Optimizer; seconds = 3600*5))
RepresentativeDaysFinders.create_plots(dft)
