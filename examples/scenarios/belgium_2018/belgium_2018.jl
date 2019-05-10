using Revise
using RepresentativeDaysFinders
using JuMP
# using GLPK
# using Gurobi
##################################################################################
# Specify location of config-file
##################################################################################
config_file = joinpath(@__DIR__,"belgium_2018.yaml")

##
using Cbc
dft = findRepresentativeDays(config_file, with_optimizer(Cbc.Optimizer; seconds = 60))

RepresentativeDaysFinders.create_plots(dft)

# Clean up data
ld = CSV.read("load_be_2018_unprocessed.csv")
# Unfortunately, this is processed as strings - Why Julia, WHY???
ldp = DataFrame(Timestep = Int[], Load = Union{Float64,Missing}[])
for i = 1:size(ld)[1]
    try
        push!(ldp, [i parse(Float64,ld[i,3])])
    catch
        push!(ldp, [i missing])
    end
end

describe(ldp, :nmissing) # Determine if anything missing, but not so all good.
for i = 1:size(ldp)[1]
    if ldp.Load[i] === missing
        try
            ldp.Load[i] = (ldp.Load[i-1] + ldp.Load[i+1])/2
        catch
            error("Your data is shit fam.")
        end
    end
end
describe(ldp, :nmissing) # Determine if anything missing, but not so all good.
