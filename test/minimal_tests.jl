# Only test essentials, used for testing different julia versions
using Test
ENV["GKSwstype"]="nul" #fix gr() problems in headless mode 

include(joinpath(@__DIR__, "test_util.jl"))
include(joinpath(@__DIR__, "test_days_finder.jl"))
include(joinpath(@__DIR__, "test_saving_and_loading.jl"))

# Remove the files
rm(joinpath(@__DIR__, "results"), recursive=true)