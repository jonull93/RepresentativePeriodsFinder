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
# interpolation  and preprocessing done here
# df_LF_windOffshore = CSV2DataFrame(joinpath(@__DIR__, "windforecast_Offshore_new.csv"); delim = ';')
# dformat = Dict(:fstr => "y-m-d HH:MM:SSz",
#                :startindex => 1,
#                :endindexOffset => 0)
# df_LF_windOffshore_filled = RepresentativeDaysFinders.check_temporal_consistency(df_LF_windOffshore, 0.25, dformat, :Column1)
# df_LF_windOffshore_interpolated = interpolatedataframe(df_LF_windOffshore_filled,:DateFormated,:LoadFactor, :LoadFactor_windOffshore, 1.)
#
# df_LF_windOnshore = CSV2DataFrame(joinpath(@__DIR__, "wind_onshore_manual.csv"); delim = ';')
# dformat = Dict(:fstr => "y-m-d HH:MM:SSz",
#                :startindex => 1,
#                :endindexOffset => 0)
# df_LF_windOnshore_filled = RepresentativeDaysFinders.check_temporal_consistency(df_LF_windOnshore, 0.25, dformat, :Column1)
# df_LF_windOnshore_interpolated = interpolatedataframe(df_LF_windOnshore_filled,:DateFormated,:LoadFactor, :LoadFactor_windOnshore, 1.)
#
# df_LF_solar = CSV2DataFrame(joinpath(@__DIR__, "solarforecast_new.csv"); delim = ';')
# dformat = Dict(:fstr => "y-m-d HH:MM:SS",
#                :startindex => 1,
#                :endindexOffset => 0)
# RepresentativeDaysFinders.scale_df_column!(df_LF_solar,:LoadFactor, 0.01) #Loadfactor seems to be given in percent
#
# df_LF_solar_filled = RepresentativeDaysFinders.check_temporal_consistency(df_LF_solar, 0.25, dformat, :Column1)
# df_LF_solar_interpolated = interpolatedataframe(df_LF_solar_filled,:DateFormated,:LoadFactor, :LoadFactor_PV, 1.)
#
#
# df_Load = CSV2DataFrame(joinpath(@__DIR__, "Total Load - Day Ahead _ Actual_2017BE.csv"); delim = ',')
# dformat = Dict(:fstr => "d.m.y HH:MM",
#                :startindex => 1,
#                :endindexOffset => 19)
#
# df_Load_filled = RepresentativeDaysFinders.check_temporal_consistency(df_Load, 0.25,dformat, Symbol("Time (CET)"))
# df_Load_interpolated = interpolatedataframe(df_Load_filled,:DateFormated,Symbol("Day-ahead Total Load Forecast [MW] - BZN|BE"), :DA_Load, 1.)
#
# LF_joint = join(df_LF_windOffshore_interpolated,df_LF_windOnshore_interpolated, on = [(:DateFormated)],makeunique=true)
# LF_joint = join(LF_joint,df_LF_solar_interpolated, on = [(:DateFormated)],makeunique=true)
# LF_joint_Final = join(LF_joint,df_Load_interpolated, on = [(:DateFormated)],makeunique=true)
#
# CSV.write(joinpath(@__DIR__, "data_interpolated.csv"), LF_joint_Final)

##################################################################################
# Specify location of config-file
##################################################################################
# using Cbc
using Gurobi
# using CPLEX
config_file = joinpath(@__DIR__,"DA_LOAD_forcast_BE_2018_load_only.yaml")
dft = findRepresentativeDays(config_file, with_optimizer(Gurobi.Optimizer; TimeLimit = 60))
# dft = findRepresentativeDays(config_file, with_optimizer(Cbc.Optimizer; seconds = 60*3*60))
# dft = findRepresentativeDays(config_file, with_optimizer(CPLEX.Optimizer; CPX_PARAM_TILIM = 60*3))
RepresentativeDaysFinders.create_plots(dft)
