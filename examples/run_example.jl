using Revise

using RepresentativeDaysFinders
using JuMP
# using GLPK
# using Gurobi
##################################################################################
# Specify location of config-file
##################################################################################
config_file = normpath(joinpath(@__DIR__, "scenarios", "test_study.yaml"))


##
using Cbc
dft = findRepresentativeDays(config_file, with_optimizer(Cbc.Optimizer; seconds = 60))
# findRepresentativeDays(config_file, with_optimizer(GLPK.Optimizer; presolve=true, msg_lev=GLPK.MSG_ALL, tm_lim=180*1000))
# Juno.@enter findRepresentativeDays(config_file, with_optimizer(GLPK.Optimizer; presolve=true, msg_lev=GLPK.MSG_ALL, tm_lim=180*1000))

RepresentativeDaysFinders.create_plots(dft)
