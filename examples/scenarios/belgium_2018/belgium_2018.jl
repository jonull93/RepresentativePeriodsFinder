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

##
using Cbc
# dft = findRepresentativeDays(config_file, with_optimizer(Cbc.Optimizer; seconds = 60))

# RepresentativeDaysFinders.create_plots(dft)
