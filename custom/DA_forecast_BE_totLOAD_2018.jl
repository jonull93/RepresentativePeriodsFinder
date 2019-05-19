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
#interpolation  and preprocessing done here
df_LF_wind = CSV2DataFrame(joinpath(@__DIR__, "windforecast_ELIA.csv"); delim = ';')
dformat = Dict(:fstr => "y-m-d HH:MM:SSz",
               :startindex => 1,
               :endindexOffset => 0)
df_LF_wind_filled = RepresentativeDaysFinders.check_temporal_consistency(df_LF_wind, 0.25, dformat, :Column1)
df_LF_wind_interpolated = interpolatedataframe(df_LF_wind_filled,:DateFormated,:LoadFactor, :LoadFactor_wind, 1.)

df_LF_solar = CSV2DataFrame(joinpath(@__DIR__, "solarforecast_ELIA.csv"); delim = ';')
dformat = Dict(:fstr => "y-m-d HH:MM:SS",
               :startindex => 1,
               :endindexOffset => 0)
df_LF_solar_filled = RepresentativeDaysFinders.check_temporal_consistency(df_LF_solar, 0.25, dformat, :Column1)
df_LF_solar_interpolated = interpolatedataframe(df_LF_solar_filled,:DateFormated,:LoadFactor, :LoadFactor_PV, 1.)


df_Load = CSV2DataFrame(joinpath(@__DIR__, "Total Load - Day Ahead _ Actual_2017BE.csv"); delim = ',')
dformat = Dict(:fstr => "d.m.y HH:MM",
               :startindex => 1,
               :endindexOffset => 19)

df_Load_filled = RepresentativeDaysFinders.check_temporal_consistency(df_Load, 0.25,dformat, Symbol("Time (CET)"))
df_Load_interpolated = interpolatedataframe(df_Load_filled,:DateFormated,Symbol("Day-ahead Total Load Forecast [MW] - BZN|BE"), :DA_Load, 1.)

LF_joint = join(df_LF_wind_interpolated,df_LF_solar_interpolated, on = [(:DateFormated)],makeunique=true)
LF_joint_Final = join(LF_joint,df_Load_interpolated, on = [(:DateFormated)],makeunique=true)

CSV.write(joinpath(@__DIR__, "data_interpolated.csv"), LF_joint_Final)

##################################################################################
# Specify location of config-file
##################################################################################
using Cbc
config_file = joinpath(@__DIR__,"DA_LOAD_forcast_BE_2018.yaml")
dft = findRepresentativeDays(config_file, with_optimizer(Cbc.Optimizer; seconds = 3600*5))
RepresentativeDaysFinders.create_plots(dft)
