using RepresentativeDaysFinders: DaysFinderTool
using Test
using YAML



@testset "RepresentativeDaysFinders" begin
    config_file = joinpath(@__DIR__, "data", "unit_testing.yaml")
    global dft = DaysFinderTool(config_file)

    @testset "TimeSeries" begin

        include("time_series/test_TimeSeries.jl")

    end


end
