# Only test essentials, used for testing different julia versions
using Test

try
    using Pkg
    Pkg.remove("HiGHS") # To fix issue in 1.2
catch e
    @warn e
end

include(joinpath(@__DIR__, "test_util.jl"))
include(joinpath(@__DIR__, "test_days_finder.jl"))
include(joinpath(@__DIR__, "test_saving_and_loading.jl"))
