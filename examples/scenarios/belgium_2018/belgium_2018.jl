using Revise
using RepresentativeDaysFinders
using JuMP
using CSV
using DataFrames
using Interpolations
using Cbc

# Electrical load data
csv_file = "load_be_2018_unprocessed.csv"
D = ENTSOEcsv2dataframe(csv_file, 3, :Load)
D = interpolatedataframe(D, :Load, 0.25, 1.0)

# Renewables data
csv_file = "generation_be_2018_unprocessed.csv"
RP_WindOnshore = ENTSOEcsv2dataframe(csv_file, 23, :OnshWindGen)
RP_WindOnshore = interpolatedataframe(RP_WindOnshore, :OnshWindGen)
RP_PV = ENTSOEcsv2dataframe(csv_file, 20, :PV)
RP_PV = interpolatedataframe(RP_PV, :PV)

# Residual load data
DRL = DataFrame(Timestep = D.Timestep, ResidualLoad = D.Load - RP_PV.PV - RP_WindOnshore.OnshWindGen)

# Normalise wind by total capacity
RP_WindOnshore = normalisegeneration(RP_WindOnshore, :OnshWindGen, 1979)
RP_PV = normalisegeneration(RP_PV, :PV, 3369)

# Combine the data frames and save to a .csv file
data = DataFrame(Timestep = D.Timestep, Load = D.Load, LFS = RP_PV.PV, LFW = RP_WindOnshore.OnshWindGen, ResidualLoad = DRL.ResidualLoad)
data = data[1:end-1,:] # Because otherwise it's 1 day too long...

# Write to .csv file
CSV.write("Belgium_2018_data.csv",data)

# Run representative days finder
config_file = joinpath(@__DIR__,"belgium_2018.yaml")
dft = findRepresentativeDays(config_file, with_optimizer(Cbc.Optimizer))

# Create plots
RepresentativeDaysFinders.create_plots(dft)
