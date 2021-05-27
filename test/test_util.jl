using RepresentativePeriodsFinder
RPF = RepresentativePeriodsFinder

# TODO: Test basic stuff like getting the number of representative periods
# to be found and stuff like that

@testset "Utilities" begin
    pf = PeriodsFinder()
    @test typeof(pf) <: RPF.PeriodsFinder

    # TODO: test has ordering error makes sense
    pf.config["method"] = Dict(
        "options" => Dict(),
        "optimization" => Dict(
            "time_series_error" => Dict(
                "weight" => 0.5,
                "type" => "squared"
            )
        )
    )
    @test RPF.has_ordering_error(pf) == true
    
    pf.config["method"]["options"]["ordering_error"] = Dict(
        "ord_err_1" => Dict(
            "function" => x -> 1,
            "weight" => 0.0
        )
    )
    @test RPF.has_ordering_error(pf) == true
    
    delete!(pf.config["method"]["optimization"], "time_series_error")
    @test RPF.has_ordering_error(pf) == true
end