using Revise
using RepresentativeDaysFinders
using JuMP
using CSV
using DataFrames
using Interpolations
# using GLPK
# using Gurobi
##################################################################################
# Specify location of config-file
##################################################################################
config_file = joinpath(@__DIR__,"belgium_2018.yaml")
csv_file = "load_be_2018_unprocessed.csv"
df = ENTSOEcsv2dataframe(csv_file, 3, :Load)
df = interpolatedataframe(df, :Load, 0.25, 1.0)

##
using Cbc
# dft = findRepresentativeDays(config_file, with_optimizer(Cbc.Optimizer; seconds = 60))

# RepresentativeDaysFinders.create_plots(dft)
